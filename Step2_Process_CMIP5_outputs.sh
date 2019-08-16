#!/bin/bash

#PBS -P dt6
#PBS -l walltime=3:00:00
#PBS -l mem=5000MB
#PBS -l ncpus=1
#PBS -j oe
#PBS -q normal
#PBS -l wd
#PBS -l other=gdata1
#PBS -l jobfs=1GB


### REQUIRES: ###
# 1) CDO command line tools
# 2) NCO
# 3) R with packages raster and ncdf4
# 4) python with package pandas


# Code processes monthly CMIP5/CMIP6 outputs
# i) merges all years into one file (when applicable),
# ii) converts required variables to desired units,
# iii) extracts desired time period
# iiii) masks out non-land grid cells (fi mask_oceans selected)
# iiiii) regrids data from Pacific- to Greenwich-central


# Note:
# Code will look for land masks in this folder:
# ${IN_DIR}/../Processed_masks/${E}/${M}/
# where E is experiment and M model
# Make sure to set "get_land_masks <- TRUE" in Step 1 if using mask_oceans=true. 
# If no mask available, model is skipped



module load cdo
module load nco
module load R

#This needs fixing on raijin, have to load miniconda to use pandas
source activate /home/561/amu561/miniconda2
#module load python



####################
### 1. SET PATHS ###
####################

#Direcotry for storing processed datasets
DIR="/g/data1/w35/amu561/CMIP6_drought/CMIP6_Data"



######################
### 2. SET OPTIONS ###
######################

dataset="cmip6"

#Clef search with options for models, experiments, variables etc.
# search_criteria="--local $dataset --experiment historical --experiment ssp585 \
#                  --experiment ssp245 --variable mrro --variable mrros \
#                  --variable pr --table Lmon --table Amon"


#for testing DELETE LATER
search_criteria="--local $dataset --experiment historical --variable mrro --table Lmon"


#Set start and end year 
year_start=1850
year_end=2100


#Mask oceans? Set to true (masking) or false (no masking). If set to true and no
#mask is found, the model and variable is skipped.
mask_oceans=true


#-------------------------------------------------------------------------------

#end of user defined settings


####################
### Create PATHS ###
####################


#Directory for saving symbolic links to original data
IN_DIR=$DIR"/CMIP6_Data"
mkdir -p $IN_DIR

#Directory for processed folder
OUT_DIR=$DIR"/Processed_CMIP6_data/"
mkdir -p $OUT_DIR

#Temporary directory for storing interim search results
TEMP_DIR=$DIR"/temp_res"
mkdir -p $TEMP_DIR



######################################
### Search database to find models ###
######################################

#Search database to find all available data
clef $search_criteria >> $TEMP_DIR/"cmip6_clef_search_results.csv"

#Filter search results to find common models and ensemble members


Rscript 



########################
### Process datasets ###
########################


#Find experiments
experiments=`ls $IN_DIR`


#Loop through experiments
for E in $experiments
do


    echo "Processing experiment: ${E}"


    #Find variables (var_table/var combination)
    vars=`find ${IN_DIR}/${E} -maxdepth 2 -mindepth 2 ! -path . -type d`


    #Loop through variables
    for V in $vars
    do


        echo "Processing variable: ${V}"

        #Variable name without path (used to create output file)
        var_short=`echo "$V" | sed 's!.*/!!'`

        #Replace nep with nee in output structure (converting NEP variable to NEE below, keeping consistent with this)
        if [ $var_short == "nep" ]; then
            var_short="nee"
        fi


        #Find models
        models=`find $V -maxdepth 1 -mindepth 1 ! -path . -type d | sed 's!.*/!!'`



        #Loop through models
        for M in $models
        do

            echo "Processing model : ${M}"

            #Find ensemble name
            ens=`ls $V/$M`

            #find files for model+ensemble
            files=`find $V/$M/$ens/*.nc`

            #Determine number of files
            no_files=`find $V/$M/$ens/*.nc | wc -l`


            #Input file to cdo functions (replaced below if files need merging)
            in_file=$files
            
            
            #If selected to mask oceans, first check if can find a land mask file
            #if not, skip model
  
            if $mask_oceans; then
              
              #Find land mask file
              mask_file=`find ${IN_DIR}/../Processed_masks/${E}/${M}/ -name "*${M}*.nc"`

              #If couldn't find mask file, skip model
              if [[ -z $mask_file ]]; then
                echo "WARNING: Could not find mask for ${M}, skipping model"
                continue
              fi

            fi
            

            #Create output path
            processed_path="${OUT_DIR}/${E}/${var_short}/${M}/"

            mkdir -p $processed_path


            #Skips files already processed to speed up processing
            ext_files=(${processed_path}/*regrid.nc)
            
            if [ -e ${ext_files[0]} ]; then
              echo "${M} - ${var_short} already processed, skipping"
              continue
            fi


            ### If multiple files, merge files ###

            if [ `echo $no_files` -gt 1 ]
            then

                #Create output file name
                out_merged="${processed_path}/Merged_${var_short}_${M}_${E}_${ens}.nc"

                #Set so ignores overlapping time periods (if present)
                export SKIP_SAME_TIME=1

                #Merge time
                cdo -L mergetime $files $out_merged

                #Find time period (brackets turn result into array)
                years=(`cdo showyear $out_merged`)

                #Write correct years to output file name, finding first and last years (use this as input to cdo)
                in_file=${out_merged%".nc"}"_"`echo ${years[0]}`"01-"`echo ${years[${#years[@]} - 1]}`"12.nc"


                #Replace merged output file with corrected file name
                mv $out_merged $in_file

            fi


            ###############################################################
            ### Convert units, mask ocean cells and extract time period ###
            ###############################################################


            # #GFDL mask files have multiple variables, need to extract "sftlt"
            # #or the division below fails.
            # if [[ $M =~ "GFDL" ]]; then
            # 
            #     #Create file name for fixed file
            #     new_mask_file=${mask_file}"_${var_short}_fixed.nc"
            # 
            #     #Select variable sftlf
            #     cdo -L selname,'sftlf' $mask_file $new_mask_file
            # 
            #     #Replace mask file with new file
            #     mask_file=$new_mask_file
            # 
            # 
            #     #Also fix $in_file, has an extra variable
            #     fixed_in_file=${in_file}"_fixed.nc"
            # 
            #     #Select variable being processed
            #     cdo -L selname,$var_short $in_file $fixed_in_file
            # 
            #     #Remove old $in_file and replace with new fixed file
            #     rm $in_file
            #     in_file=$fixed_in_file
            # fi
            # 
            
            
            #Create output file name (complete when running each cdo command)
            processed_file="${processed_path}/${var_short}_${M}_${E}_${ens}"


            #########################
            ### Process variables ###
            #########################
            
            ###--- Mask oceans ---###
            
            
            #First mask for oceans if this option was selected
            if $mask_oceans; then
              
              #Mask ocean cells
              cdo -L div $in_file -gec,99 $mask_file  ${processed_file}_temp.nc
              
              #Replace input file for next step
              in_file=${processed_file}_temp.nc
              
            fi
              
              
          ###--- Select years ---###
          
          #Check that wanted years exist in files 
          years=(`cdo showyear $in_file`)




          #find first and last years 
          `echo ${years[0]}`
          
          `echo ${years[${#years[@]} - 1]}

          
          #If time period in file shorter, stop 
            
          
          #
          #Select years
          cdo selyear,$year_start/$year_end $in_file $out_file
  
          infile=
          
          
            

            ###--- Air temperature ---###
            if [[ $V =~ "Amon/tas" ]]; then

                #Create output file
                out_file="${processed_file}_deg_C_${year_start}_${year_end}_${E}.nc"

                #Convert Kelvin to Celsius
                cdo expr,'tas=tas-273.15' -setunit,'degrees C' $in_file $out_file


            ###--- Evapotranspiration ---###
            elif [[ $V =~ "Amon/evspsbl" ]]; then

                #Create output file
                out_file="${processed_file}_mm_month_${year_start}_${year_end}_${E}.nc"

                #Convert from kg m-2 s-1 to mm/month
                cdo -L muldpm -expr,'evspsbl=evspsbl*60*60*24' -setunit,'mm/month' $in_file $out_file


            ###--- Precipitation ---###
            elif [[ $V =~ "Amon/pr" ]]; then

                #Create output file
                out_file="${processed_file}_mm_month_${year_start}_${year_end}_${E}.nc"

                #Convert from kg m-2 s-1 to mm/month
                cdo muldpm -expr,'pr=pr*60*60*24' -setunit,'mm/month' $in_file ${processed_file}_temp1.nc

                #Replace negative rainfall with zero
                ncap2 -s 'where(pr < 0) pr=0' -O ${processed_file}_temp1.nc $out_file

                rm ${processed_file}_temp1.nc


            ###--- Gross primary production ---###
            elif [[ $V =~ "Lmon/gpp" ]]; then

                #Create output file
                out_file="${processed_file}_gC_m2_month_${year_start}_${year_end}_${E}.nc"
    
                #Convert from kg m-2 s-1 to g C m-2 month-1
                cdo muldpm -expr,'gpp=gpp*60*60*24*1000' -setunit,'g C m-2 month-1' $in_file ${processed_file}_temp1.nc

                #Replace negative GPP with zero (-O switch overwrites existing file if any)
                ncap2 -s 'where(gpp < 0) gpp=0' -O ${processed_file}_temp1.nc $out_file

                rm ${processed_file}_temp1.nc


            ###--- Net ecosystem exchange ---###
            elif [[ $V =~ "Lmon/nep" ]]; then


                #Create output file
                out_file="${processed_file}_gC_m2_month_${year_start}_${year_end}_${E}.nc"

                #Convert from kg m-2 s-1 to g C m-2 month-1, and change sign (to go from NEP to NEE)
                cdo muldpm -expr,'nep=nep*60*60*24*1000*(-1)' -setunit,'g C m-2 month-1' $in_file ${processed_file}_temp1.nc

                #Change variable name from NEP to NEE
                cdo chname,nep,nee ${processed_file}_temp1.nc $out_file

                rm ${processed_file}_temp1.nc

  

            ###--- Surface runoff ---###
            #Important: have mrros before mrro or code will use mrro (because of if argument V "contains" mrro)
            elif [[ $V =~ "Lmon/mrros" ]]; then

                #Create output file
                out_file="${processed_file}_mm_month_${year_start}_${year_end}_${E}.nc"

                #Convert from kg m-2 s-1 to mm/month
                cdo muldpm -expr,'mrros=mrros*60*60*24' -setunit,'mm/month' $in_file ${processed_file}_temp1.nc

                #Replace negative runoff with zero
                ncap2 -s 'where(mrros < 0) mrros=0' -O ${processed_file}_temp1.nc $out_file

                rm ${processed_file}_temp1.nc


            ###--- Total runoff ---###
            elif [[ $V =~ "Lmon/mrro" ]]; then

                #Create output file
                out_file="${processed_file}_mm_month_${year_start}_${year_end}_${E}.nc"

                #Convert from kg m-2 s-1 to mm/month
                cdo muldpm -expr,'mrro=mrro*60*60*24' -setunit,'mm/month' $in_file ${processed_file}_temp1.nc

                #Replace negative runoff with zero
                ncap2 -s 'where(mrro < 0) mrro=0' -O ${processed_file}_temp1.nc $out_file

                rm ${processed_file}_temp1.nc



            ###--- Daily tasmax ---###
            elif [[ $V =~ "day/tasmax" ]]; then

                #Create output file
                out_file="${processed_file}_deg_C_${year_start}_${year_end}_${E}.nc"

                #Convert Kelvin to Celsius (change unit, average to monthly, select years and convert to C)
                cdo expr,'tasmax=tasmax-273.15' -monmean -setunit,'degrees C' $in_file $out_file


            ###--- Daily precip ---###
            elif [[ $V =~ "day/pr" ]]; then

                #Create output file
                out_file="${processed_file}_mm_day_${year_start}_${year_end}_${E}.nc"

                #Convert Kelvin to Celsius (change unit, average to monthly, select years and convert to C)
                cdo expr,'pr=pr*60*60*24' -setunit,'mm/day' $in_file $out_file                
                

            ###--- All other variables ---###
            else

                #Create output file
                out_file="${processed_file}_${year_start}_${year_end}_${E}.nc"

                mv $in_file $out_file
                

            fi


            ### Tidy up ###
            
            #Remove temp files
            temp_files=`find ${processed_path} -name *temp.nc`
            if [[ -e $temp_files ]]; then rm $temp_files; fi
            
            #Remove merged files
            merged_files=`find ${processed_path} -name Merged_*.nc`
            if [[ -e $merged_files ]]; then rm $merged_files; fi

            #Remove fixed GFDL mask file created earlier
            if [[ $M =~ "GFDL" ]]; then rm $new_mask_file; fi
            
              
                              
            ###################################################################
            ###--- Regrid output file from Pacific- to Greenwich-central ---###
            ###################################################################

            # Modify output file name
            outfile_regrid=${out_file%".nc"}"_regrid.nc"

		        #Regrid using CDO sellonlatbox (note this adjusts lon variable
            #to range [-180, 180] but not lon_bounds)
		        cdo sellonlatbox,-180,180,-90,90 $out_file $outfile_regrid

            #Above CDO command doesn't fix lon_bounds, using a python
            #script to fix these (provided by Arden)
            python fix_lon.py "${outfile_regrid}"

		        #Remove temp file if cropping
		        #rm $outfile_temp

            

			      #Remove non-regridded file and merged time file
			      rm $out_file



            ###########################################
            ###--- Data quality and error checks ---###
            ###########################################

            #Sanity check, does output file exist?
            files_existing=`ls $processed_path`

            #Check if empty string. If so, cat to file
            if [ -z "$files_existing" ]; then

               echo $processed_path >> failed_files.txt

            fi


            #Calculate and save variable global mean. Use to check that variables are not corrupted
            check_dir="${OUT_DIR}/Data_checks/${E}/${var_short}/Global_mean/"
            mkdir -p $check_dir
            
            #The python script (fix_lon.py) resaves the file with a different name (setgrid*)
            #Use this when calculating global mean
            outfile_setgrid=${out_file%".nc"}"_regrid_setgrid.nc"
            
			      cdo fldmean $outfile_regrid ${check_dir}/${M}_global_mean_${var_short}.nc


            
            ##############################
            ###--- Data check cont. ---###
            ##############################

            #Plot monthly mean of variable in all models (for regridded files)
            #Used to check no errors have occurred during processing (or that original file is not corrupted)

            plot_dir="${OUT_DIR}/Data_checks/Plots/${M}/"
            mkdir -p $plot_dir

            #Create R script for plotting
            cat > R_plot.R << EOF
            library(raster)
            library(ncdf4)


            ### Map of first time slice ###
            files_regrid <- list.files(path="${OUT_DIR}/${E}/${var_short}/${M}/", recursive=TRUE, 
                                       pattern="setgrid.nc", full.names=TRUE)    #regridded
                                       
            data_regrid <- lapply(files_regrid, brick, stopIfNotEqualSpaced=FALSE)


            pdf("${plot_dir}/${var_short}_${E}_${M}_monthly_mean_regridded.pdf", 
                height=5, width=8)
            par(mai=c(0.2,0.2,0.2,0.6))
            par(mfcol=c(ceiling(sqrt(length(data_regrid))), ceiling(sqrt(length(data_regrid)))))
            lapply(data_regrid, function(x) plot(mean(x), ylab="", xlab="", yaxt="n", xaxt="n"))
            dev.off()


        		### Global mean time series ###
        		files_mean <- list.files(path="${check_dir}", recursive=TRUE, 
                                     pattern="${M}_global_mean", full.names=TRUE)

        		nc_handles <- lapply(files_mean, nc_open)
        		nc_data    <- lapply(nc_handles, ncvar_get, varid="${var_short}")

        		pdf("${plot_dir}/${var_short}_${E}_${M}_global_mean_timeseries.pdf", 
                height=13, width=40)
        		par(mai=c(0.2,0.2,0.2,0.6))
        		par(mfcol=c(ceiling(sqrt(length(nc_data))), ceiling(sqrt(length(nc_data)))))
        		lapply(nc_data, function(x) plot(x, type="l", col="blue", ylab="", xlab="", 
                   yaxt="n", xaxt="n"))
        		dev.off()


EOF

            #Run script and remove afterwards
            Rscript R_plot.R
            rm R_plot.R


        done #models

    done #vars

done #experiments
























