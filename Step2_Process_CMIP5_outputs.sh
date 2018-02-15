#!/bin/bash

#PBS -P dt6
#PBS -l walltime=17200
#PBS -l mem=5000MB
#PBS -l ncpus=1
#PBS -j oe
#PBS -q normal
#PBS -l wd
#PBS -l other=gdata1
#PBS -l jobfs=1GB


# Code processes monthly CMIP5 outputs
# i) merges all years into one file (when applicable),
# ii) converts required variables to desired units,
# iii) extracts desired time period
# iiii) masks out non-land grid cells
# iiiii) regrids data from Pacific- to Greenwich-central (note longitude still
#        ranges 0-360 as not sure how to fix this in NCL, use Step3 R code to fix this)


# NB.    !!!!!!!!!!!!!!!!!
# Code will look for land masks in this folder:
# ${IN_DIR}/../Processed_masks/${E}/${M}/
# where E is experiment and M model
# Make sure to save land masks here first


# WARNING !!!!! : manually handles CMCC-CESM gpp and nep data, setting a 
#                 symbolic link to a manually corrected file


module load cdo
module load nco
module load ncl
module load R


####################
### 1. SET PATHS ###
####################

#Folder used in Step 1 
IN_DIR="/g/data1/w35/amu561/CMIP5_fluxnet/Data"

#Desired output folder
OUT_DIR="/g/data1/w35/amu561/CMIP5_fluxnet/Processed_data_1950_2099/"


######################
### 2. SET OPTIONS ###
######################

#Set start and end year 
year_start=1950
year_end=2099



#-------------------------------------------------------------------------------

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


        echo "ECHOING models : ${models}"

        #Loop through models
        for M in $models
        do


            #Find ensemble name
            ens=`ls $V/$M`

            #find files for model+ensemble
            files=`find $V/$M/$ens/*.nc`

            #Determine number of files
            no_files=`find $V/$M/$ens/*.nc | wc -l`


            #Input file to cdo functions (replaced below if files need merging)
            in_file=$files


            #Create output path
            processed_path="${OUT_DIR}/${continent}/${var_short}/${M}/"

            mkdir -p $processed_path



            ### If multiple files, merge files ###

            if [ `echo $no_files` -gt 1 ]
            then

                #Create output file name
                out_merged="${processed_path}/Merged_${var_short}_${M}_${E}_${ens}.nc"

                #Set so ignores overlapping time periods (if present)
                export SKIP_SAME_TIME=1

                #Merge time
                cdo mergetime $files $out_merged

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


            #Find land mask file
            mask_file=`find ${IN_DIR}/../Processed_masks/${E}/${M}/ -name "*.nc"`


            #GFDL files have multiple variables, need to extract "sftlt"
            #or the division below fails.
            if [[ $M =~ "GFDL" ]]; then
              
                #Create file name (use continent and varname so can run the code
                #multiple times at the same time without inteference)
                new_mask_file=${mask_file}"_${continent}_${var_short}_fixed.nc"
                
                #Select variable sftlf
                cdo selname,'sftlf' $mask_file $new_mask_file
                
                #Replace mask file with new file
                mask_file=$new_mask_file
                
                
                #Also fix $in_file, has an extra variable
                fixed_in_file=${in_file}"_fixed.nc"
                
                #Select variable being processed
                cdo selname,$var_short $in_file $fixed_in_file
                
                #Remove old $in_file and replace with new fixed file
                rm $in_file
                in_file=$fixed_in_file
            fi
            
            
            #Create output file name (complete when running each cdo command)
            processed_file="${processed_path}/${var_short}_${M}_${E}_${ens}"


            echo "ECHOING processed file : ${processed_file}"



            ### Process variables ###


            ###--- Air temperature ---###
            if [[ $V =~ "Amon/tas" ]]; then

                #Create output file
                out_file="${processed_file}_deg_C_${year_start}_${year_end}_${continent}.nc"

                #Mask ocean cells
                cdo div $in_file -gec,99 $mask_file  ${processed_file}_temp.nc

                #Convert Kelvin to Celsius
                cdo expr,'tas=tas-273.15' -selyear,$year_start/$year_end -setunit,'degrees C' ${processed_file}_temp.nc $out_file

                rm ${processed_file}_temp.nc



            ###--- Evapotranspiration ---###
            elif [[ $V =~ "Amon/evspsbl" ]]; then


                #Create output file
                out_file="${processed_file}_mm_month_${year_start}_${year_end}_${continent}.nc"

                #Mask ocean cells
                cdo div $in_file -gec,99 $mask_file  ${processed_file}_temp.nc


                #These models have ET 1000 times too high, divide ET to fix (see CMIP errata, http://cmip-pcmdi.llnl.gov/cmip5/errata/cmip5errata.html)
                if [[ $M =~ "NorESM1-M" ]]; then
                    #Convert from kg m-2 s-1 to mm/month, and divide by 1000 to fix corrupted file
                    cdo muldpm -expr,'evspsbl=evspsbl*60*60*24/1000' -selyear,$year_start/$year_end -setunit,'mm/month' ${processed_file}_temp.nc $out_file

                #This model has wrong ET sign (negative ET), change sign by multiplying with -1
                elif [[ $M =~ "CMCC-C" ]]; then
                    cdo muldpm -expr,'evspsbl=evspsbl*60*60*24*(-1)' -selyear,$year_start/$year_end -setunit,'mm/month' ${processed_file}_temp.nc $out_file

                else
                    #Convert from kg m-2 s-1 to mm/month
                    cdo muldpm -expr,'evspsbl=evspsbl*60*60*24' -selyear,$year_start/$year_end -setunit,'mm/month' ${processed_file}_temp.nc $out_file
                fi



                rm ${processed_file}_temp.nc



            ###--- Precipitation ---###
            elif [[ $V =~ "Amon/pr" ]]; then

                #Create output file
                out_file="${processed_file}_mm_month_${year_start}_${year_end}_${continent}.nc"

                #Mask ocean cells
                cdo div $in_file -gec,99 $mask_file  ${processed_file}_temp.nc

                #Convert from kg m-2 s-1 to mm/month
                cdo muldpm -expr,'pr=pr*60*60*24' -selyear,$year_start/$year_end -setunit,'mm/month' ${processed_file}_temp.nc ${processed_file}_temp1.nc

                #Replace negative rainfall with zero
                ncap2 -s 'where(pr < 0) pr=0' -O ${processed_file}_temp1.nc $out_file


                rm ${processed_file}_temp.nc
                rm ${processed_file}_temp1.nc



            ###--- Gross primary production ---###
            elif [[ $V =~ "Lmon/gpp" ]]; then

                #Create output file
                out_file="${processed_file}_gC_m2_month_${year_start}_${year_end}_${continent}.nc"


                #If CMCC model, link to a manually corrected file (output corrupted)
                if [[ $M =~ "CMCC-CESM" ]]; then

                    rm $out_file #remove in case exists, copy will fail otherwise
                    cp -s ${IN_DIR}/../Files_to_replace_corrupted/gpp_Lmon_CMCC-CESM_historical_r1i1p1_1901-2004_v20140417_manually_fixed_monthly_total.nc ${processed_file}_temp.nc

                    #Select appropriate years
                    cdo selyear,$year_start/$year_end ${processed_file}_temp.nc $out_file

                    rm ${processed_file}_temp.nc

                else

                    #Mask ocean cells
                    cdo div $in_file -gec,99 $mask_file  ${processed_file}_temp.nc

                    #Convert from kg m-2 s-1 to g m-2 s-1
                    cdo muldpm -expr,'gpp=gpp*60*60*24*1000' -selyear,$year_start/$year_end -setunit,'g C m-2 month-1' ${processed_file}_temp.nc $  {processed_file}_temp1.nc

                    #Replace negative rainfall with zero (-O switch overwrites existing file if any)
                    ncap2 -s 'where(gpp < 0) gpp=0' -O ${processed_file}_temp1.nc $out_file


                    rm ${processed_file}_temp.nc
                    rm ${processed_file}_temp1.nc

                fi


            ###--- Net ecosystem exchange ---###
            elif [[ $V =~ "Lmon/nep" ]]; then


                #Create output file
                out_file="${processed_file}_gC_m2_month_${year_start}_${year_end}_${continent}.nc"


                #If CMCC model, link to a manually corrected file (output corrupted)
                if [[ $M =~ "CMCC-CESM" ]]; then

                    rm $out_file #remove in case exists, copy will fail otherwise
                    cp -s ${IN_DIR}/../Files_to_replace_corrupted/nee_Lmon_CMCC-CESM_historical_r1i1p1_1901-2004_v20140417_manually_fixed_monthly_total.nc ${processed_file}_temp.nc


                    #Select appropriate years
                    cdo selyear,$year_start/$year_end ${processed_file}_temp.nc $out_file

                    rm ${processed_file}_temp.nc


                else

                    #Mask ocean cells
                    cdo div $in_file -gec,99 $mask_file  ${processed_file}_temp.nc

                    #Convert from kg m-2 s-1 to g m-2 s-1, and change sign (to go from NEP to NEE)
                    cdo muldpm -expr,'nep=nep*60*60*24*1000*(-1)' -selyear,$year_start/$year_end -setunit,'g C m-2 month-1' ${processed_file}_temp.nc ${processed_file}_temp1.nc

                    #Change variable name from NEP to NEE
                    cdo chname,nep,nee ${processed_file}_temp1.nc $out_file

                    rm ${processed_file}_temp.nc
                    rm ${processed_file}_temp1.nc

                fi


            ###--- Surface runoff ---###
            #Important: have mrros before mrro or code will use mrro (because of if argument V "contains" mrro)
            elif [[ $V =~ "Lmon/mrros" ]]; then

                #Create output file
                out_file="${processed_file}_mm_month_${year_start}_${year_end}_${continent}.nc"

                cdo div $in_file -gec,99 $mask_file  ${processed_file}_temp.nc

                #Convert from kg m-2 s-1 to mm/month
                cdo muldpm -expr,'mrros=mrros*60*60*24' -selyear,$year_start/$year_end -setunit,'mm/month' ${processed_file}_temp.nc ${processed_file}_temp1.nc

                #Replace negative runoff with zero
                ncap2 -s 'where(mrros < 0) mrros=0' -O ${processed_file}_temp1.nc $out_file

                rm ${processed_file}_temp.nc
                rm ${processed_file}_temp1.nc


            ###--- Total runoff ---###
            elif [[ $V =~ "Lmon/mrro" ]]; then

                #Create output file
                out_file="${processed_file}_mm_month_${year_start}_${year_end}_${continent}.nc"

                cdo div $in_file -gec,99 $mask_file  ${processed_file}_temp.nc

                #Convert from kg m-2 s-1 to mm/month
                cdo muldpm -expr,'mrro=mrro*60*60*24' -selyear,$year_start/$year_end -setunit,'mm/month' ${processed_file}_temp.nc ${processed_file}_temp1.nc

                #Replace negative runoff with zero
                ncap2 -s 'where(mrro < 0) mrro=0' -O ${processed_file}_temp1.nc $out_file

                rm ${processed_file}_temp.nc
                rm ${processed_file}_temp1.nc



            ###--- Tasmax ---###
            elif [[ $V =~ "day/tasmax" ]]; then

                #Create output file
                out_file="${processed_file}_deg_C_${year_start}_${year_end}_${continent}.nc"

                #Mask ocean cells
                cdo div $in_file -gec,99 $mask_file  ${processed_file}_temp.nc

                #Convert Kelvin to Celsius (change unit, average to monthly, select years and convert to C)
                cdo expr,'tasmax=tasmax-273.15' -selyear,$year_start/$year_end -monmean -setunit,'degrees C' ${processed_file}_temp.nc $out_file

                rm ${processed_file}_temp.nc



            ###--- Tasmax ---###
            elif [[ $V =~ "day/pr" ]]; then

                #Create output file
                out_file="${processed_file}_mm_day_${year_start}_${year_end}_${continent}.nc"

                #Mask ocean cells
                cdo div $in_file -gec,99 $mask_file  ${processed_file}_temp.nc

                #Convert Kelvin to Celsius (change unit, average to monthly, select years and convert to C)
                cdo expr,'pr=pr*60*60*24' -selyear,$year_start/$year_end -setunit,'mm/day' ${processed_file}_temp.nc $out_file

                rm ${processed_file}_temp.nc
                
                

            ###--- All other variables ---###
            else

                out_file="${processed_file}_${year_start}_${year_end}_${continent}.nc"

                cdo div $in_file -gec,99 $mask_file  ${processed_file}_temp.nc

                cdo selyear,$year_start/$year_end ${processed_file}_temp.nc $out_file

                rm ${processed_file}_temp.nc

            fi


            #Remove fixed GFDL mask file created earlier
            if [[ $M =~ "GFDL" ]]; then  
                rm $new_mask_file
            fi
            
                
            ############################################################################
            ###--- Crop and regrid output file from Pacific- to Greenwich-central ---###
            ############################################################################

            # Modify output file name
            outfile_regrid=${out_file%".nc"}"_regrid.nc"


            #If CMCC model, link to a manually corrected file (output corrupted)
            if [[ $M =~ "CMCC-CESM" && $V =~ "Lmon/gpp" ]]
            then


                rm $outfile_regrid #remove in case exists, copy will fail otherwise
                cp -s ${IN_DIR}/../Files_to_replace_corrupted/gpp_Lmon_CMCC-CESM_historical_r1i1p1_1901-2004_v20140417_manually_fixed_monthly_total_regrid.nc ${out_file%".nc"}"_temp.nc"

                #Select years
                cdo sellonlatbox,$ext -selyear,$year_start/$year_end ${out_file%".nc"}"_temp.nc" $outfile_regrid

                rm ${out_file%".nc"}"_temp.nc"


            elif [[ $M =~ "CMCC-CESM" && $V =~ "Lmon/nep" ]]
            then

                rm $outfile_regrid #remove in case exists, copy will fail otherwise
                cp -s ${IN_DIR}/../Files_to_replace_corrupted/nee_Lmon_CMCC-CESM_historical_r1i1p1_1901-2004_v20140417_manually_fixed_monthly_total_regrid.nc {out_file%".nc"}"_temp.nc"

                #Select years
                cdo sellonlatbox,$ext -selyear,$year_start/$year_end ${out_file%".nc"}"_temp.nc" $outfile_regrid

                rm ${out_file%".nc"}"_temp.nc"


            #All other models
            else

				        #Temp file name
				        outfile_temp=${out_file%".nc"}"_temp.nc"

				        #Set inputs
                export INPUTFILE=$out_file
                export OUTPUTFILE=$outfile_temp

                #Run NCL code to regrid (cd to directory where NCL codes are.
                #Not ideal but can't find linked ncl function otherwise, not sure how to fix.
                cd ${IN_DIR}"/../scripts/Bash/functions/"
                ncl fix_cmip_lon.ncl


				        #Then crop to extent
				        cdo sellonlatbox,$ext $outfile_temp $outfile_regrid


				        #Remove temp file
				        rm $outfile_temp

            fi


			      #Remove non-regridded file and merged time file
			      rm $out_file
			      rm $in_file

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
            #(mean not properly calculated, doesn't weight cells by area. Only for guidance)
            check_dir="${OUT_DIR}/Data_checks/${E}/${var_short}/Global_mean/"
            mkdir -p $check_dir
			      cdo fldmean $outfile_regrid ${check_dir}/${M}_global_mean_${continent}.nc



        done #models



        ##############################
        ###--- Data check cont. ---###
        ##############################

        #Plot monthly mean of variable in all models (for regridded files)
        #Used to check no errors have occurred during processing (or that original file is not corrupted)

        plot_dir="${OUT_DIR}/Data_checks/Plots/"
        mkdir -p $plot_dir

        #Create R script for plotting
        cat > R_plot.R << EOF
        library(raster)
        library(ncdf4)


        ### Map of first time slice ###
        files_regrid <- list.files(path="${OUT_DIR}/${E}/${var_short}", recursive=TRUE, pattern="regrid.nc", full.names=TRUE)    #regridded
		    #files        <- list.files(path="${OUT_DIR}/${E}/${var_short}", recursive=TRUE, pattern="${end_yr}.nc", full.names=TRUE) #non-regridded

        data_regrid <- lapply(files_regrid, brick, stopIfNotEqualSpaced=FALSE)
		    #data        <- lapply(files, brick, stopIfNotEqualSpaced=FALSE)

        pdf("${plot_dir}/${var_short}_${E}_all_models_monthly_mean_regridded_${continent}.pdf", height=15, width=25)
        par(mai=c(0.2,0.2,0.2,0.6))
        par(mfcol=c(ceiling(sqrt(length(data_regrid))), ceiling(sqrt(length(data_regrid)))))
        lapply(data_regrid, function(x) plot(mean(x), ylab="", xlab="", yaxt="n", xaxt="n"))
        dev.off()

		#pdf("${plot_dir}/${var_short}_${E}_all_models_monthly_mean_non-regridded_${continent}.pdf", height=15, width=25)
		#par(mai=c(0.2,0.2,0.2,0.6))
		#par(mfcol=c(ceiling(sqrt(length(data))), ceiling(sqrt(length(data)))))
		#lapply(data, function(x) plot(mean(x), ylab="", xlab="", yaxt="n", xaxt="n"))
		#dev.off()


		### Global mean time series ###
		files_mean <- list.files(path=$check_dir, recursive=TRUE, pattern="global_mean", full.names=TRUE)

		nc_handles <- lapply(files_mean, nc_open)
		nc_data <- lapply(nc_handles, ncvar_get, varid=$var_short)

		pdf("${plot_dir}/${var_short}_${E}_all_models_global_mean_timeseries_${continent}.pdf", height=13, width=40)
		par(mai=c(0.2,0.2,0.2,0.6))
		par(mfcol=c(ceiling(sqrt(length(nc_data))), ceiling(sqrt(length(nc_data)))))
		lapply(nc_data, function(x) plot(x, type="l", col="blue", ylab="", xlab="", yaxt="n", xaxt="n"))
		dev.off()


EOF


        #Run script and remove afterwards
        Rscript R_plot.R
        rm R_plot.R


    done #vars

done #experiments
























