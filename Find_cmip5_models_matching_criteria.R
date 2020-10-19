### Find and sym-link copy CMIP5 runs to match criteria ###


#clear R environment
rm(list=ls(all=TRUE))

args = commandArgs(trailingOnly=TRUE)

####################
### Set criteria ###
####################

### 1. SET PATH ###

#Directory where to copy data
outdir <- args[1]


### 2. DECIDE FILE STRUCTURE ###

#Should model files for all experiments be
#saved in same folder (specify name of folder)?
#If want to e.g. combine historical and RCP8.5 runs,
#use this option, else set to FALSE
combine  <- as.logical(args[2]) 
dir_name <- args[3] 


### 3. DECIDE IF WANT LAND MASKS ###

#Retrieves land masks for selected models 
get_land_masks <- as.logical(args[4]) 

#name of land mask variable
mask_var  <- args[5]


#------------------------------------------------------------------------------


##########################
### Get query dataset ####
##########################


#Get Clef search results
results <- as.vector(read.csv(args[6], header=FALSE, 
			      colClasses="character"))



#### NEED TO CHANGE ONCE CLEF CHANGED !!!!!!!!!!!!!!!!!!!

#Currently having to retrieve info from paths, and hardcode indices

#Create a table for search results
sorted_results <- as.data.frame(matrix(data=NA, nrow=nrow(results), ncol=8))
colnames(sorted_results)  <- c("model", "experiment", "ensemble", "variable", 
                               "grid", "version", "time_resolution", "path")


#Separate search results
all_res <- lapply(results, function(x) strsplit(x, "/"))[[1]]


#Get models
sorted_results$model      <- sapply(all_res, function(x) x[9])

#Get experiments
sorted_results$experiment <- sapply(all_res, function(x) x[10])

#Get ensemble member
sorted_results$ensemble   <- sapply(all_res, function(x) x[14])

#Get variable
sorted_results$variable   <- sapply(all_res, function(x) x[16])


missing_versions <- NA
#Some variable names include an underscore
if (any (grepl("_", sorted_results$variable))) {
  
  vars_temp <- lapply(sorted_results$variable, function(x) strsplit(x, "_"))
  
  #Grab first results (this is the variable)
  sorted_results$variable <- sapply(vars_temp, function(x) x[[1]][1])
  
  #Some of these results also include the version, need to grab that
  missing_versions <- sapply(vars_temp, function(x) x[[1]][2])
  
  miss_ind <- which(!is.na(missing_versions))
  
}


#Get version
sorted_results$version    <- sapply(all_res, function(x) x[15])

#Replace versions with those obtained from var names above
if (!(is.na(missing_versions[1]))) sorted_results$version[miss_ind] <- missing_versions[miss_ind]


#Get time resolution
sorted_results$time_resolution <- sapply(all_res, function(x) x[13])

#Get paths
sorted_results$path <- results[,1]



#Split mask variable from results

if (get_land_masks) {
  
  mask_ind <- which(sorted_results$variable == mask_var)
  
  #Get mask results
  mask_results <- sorted_results[mask_ind,]
  
  #And the rest
  sorted_results <- sorted_results[-mask_ind,]
  
}


#FIO-ESM is provided as lower and uppercase, fix this
if (any(sorted_results$model %in% c("FIO-ESM", "fio-esm"))) {
  ind <- which(sorted_results$model %in% c("FIO-ESM", "fio-esm"))
  sorted_results$model[ind] <- "FIO-ESM"
}




##########################
### Find common models ###
##########################

#Separate results by experiment
experiments <- unique(sorted_results$experiment)

res_experiment <- lapply(experiments, function(x) sorted_results[which(sorted_results$experiment==x),])


#Find models for each experiment and variable
variables <- unique(sorted_results$variable)

models <- lapply(res_experiment, function(x) lapply(variables, function(y) x$model[which(x$variable==y)]))


#Find all available models
all_models <- unique(unlist(models))


#Find common models
common_models <- vector()
for (k in 1:length(all_models)) {
  
  #Find if model exist for all experiments of a variable
  exist <- sapply(models, function(x) all(sapply(x, function(y) any(y==all_models[k]))))
  
  
  #Save model name if exists for all variables
  if(all(exist)) common_models <- append(common_models, all_models[k])
  
}


#############################
### Find ensemble members ###
#############################

ensemble <- lapply(res_experiment, function(x) lapply(variables, function(y) x$ensemble[which(x$variable==y)]))


selected_ens <- list()
for (k in 1:length(common_models)) {
  
  #Find all ensemble members for each experiment and variable
  avail_ens <- lapply(1:length(ensemble), function(exp) mapply(function(ens,mod) ens[which(mod==common_models[k])], 
                                                               ens=ensemble[[exp]], mod=models[[exp]], SIMPLIFY=FALSE))
  
  
  #Find all unique ensembles
  all_ens <- unique(unlist(avail_ens))
  
  #Find common ensemble members
  common_ens <- vector()  
  for (e in 1:length(all_ens)) {
    
    #Find if model exist for all experiments of a variable
    exist <- sapply(avail_ens, function(x) all(sapply(x, function(y) any(y==all_ens[e]))))
    
    #Save model name if exists for all variables
    if (all(exist)) common_ens <- append(common_ens, all_ens[e])
  }
  
  # #Select first available ensemble
  # if (length(common_ens) > 0) {
  #   #Sorting this way so r1i1p1 comes before r10i1p1 etc.
  selected_ens[[k]] <- common_ens
  # }
}


# Check if didn't find any common ensemble members for a model
# Remove model in this case
if (any(sapply(selected_ens, length) < 1)) {
  rm_ind        <- which(sapply(selected_ens, length) < 1) #model(s) to remove
  common_models <- common_models[-rm_ind]
  selected_ens <- selected_ens[-rm_ind]
}




####################################
### Collate final search results ###
####################################

final_results <- sorted_results


### Remove extra models ###

mod_ind <- which(!(final_results$model %in% common_models))

if (length(mod_ind) >0) final_results <- final_results[-mod_ind,]


### Remove extra ensembles ###

for (m in 1:length(common_models)) {
  
  #Find model 
  ind <- which(final_results$model == common_models[m])
  
  #Find extra ensemble members
  ens_ind <- which(!(final_results[ind,]$ensemble %in% selected_ens[[m]]))
  
  if (length(ens_ind) >0) final_results <- final_results[-ind[ens_ind],]
  
}


### Remove extra versions ###

#Some variables might have multiple versions available,
#only use the latest

#Variables to match
match <- final_results[,c("model", "experiment", "ensemble", "variable")]

#Find duplicates
vr_ind <- which(duplicated(match))

#If found instances
if (length(vr_ind) > 0) {
  
  rm_ind <- vector()
  
  for (v in 1:length(vr_ind)) {
    
    #Finds rows that are duplicates
    common_ind <- which(apply(match, MARGIN=1, function(x) all(x == match[vr_ind[v],])))
    
    #Check which one of these is the newest version
    #(remove "v" from start, convert to numeric and find biggest)
    
    #Some models have the latest versions as "latest", not a time stamp
    #Take this is available
    if (any(final_results$version[common_ind] == "latest")) {
      
      #Take first available "latest" if several
      max_ind <- which(final_results$version[common_ind] == "latest")[1]
    
    #Else calculate from time/version stamp
    } else {
      max_ind <- which.max(as.numeric(substr(final_results$version[common_ind], 2, 9)))
      
    }
      
    
    rm_ind <- append(rm_ind, common_ind[-max_ind])
    
  }
  
  #Remove all old versions
  final_results <- final_results[-rm_ind,]
  
}




###########################################
### Create database of selected outputs ###
###########################################


#Loop through models
for (k in 1:nrow(final_results)) {
  
  #Extract data for this iteration
  entry <- final_results[k,]
  
  
  
  #If saving to same directory
  if (combine) {
    #Create output directory
    target_dir <- paste(outdir, dir_name, entry$time_resolution,
                        entry$variable, entry$model,
                        entry$ensemble, sep="/")
    
    #Else
  } else {
    #Create output directory
    target_dir <- paste(outdir, entry$experiment, entry$time_resolution,
                        entry$variable, entry$model,
                        entry$ensemble, sep="/")
  }
  
  dir.create(target_dir, recursive=TRUE, showWarnings = FALSE)
  
  
  
  
  #Basic sanity checks:
  #Check that path contains correct model, variable and experiment  
  
  #Other models
  
  #Model check
  if (!grepl(entry$model, entry$path)) {
    stop("File path doesn't contain correct model")
  } 
  #Variable check
  if (!grepl(entry$variable, entry$path)) {
    stop("File path doesn't contain correct variable")
  }
  
  
  #Experiment check
  if (!grepl(entry$experiment, entry$path)) {
    stop("File path doesn't contain correct experiment")
  }
  
  
  
  #Create symbolic link to file (need to add pattern because some models 
  #will otherwise return surplus variables. Also adding "_" to avoid similar 
  #var names to be returned)
  
  file.symlink(from=list.files(entry$path, full.names=TRUE, 
                               pattern=paste0(entry$variable, "_")), 
               to=target_dir)
  
}




#########################################
### Retrieve corresponding land masks ###
#########################################

if (get_land_masks) {
  
  for (k in 1:length(common_models)) {
    
    
    #Extract results for model
    mod_res <- mask_results[which(mask_results$model == common_models[k]), ]
    
    
    #If found results, copy file
    if (nrow(mod_res) > 0) {
      
      
      ### Check for duplicate versions ###
      
      #Check if any experiment/ensemble/variable combo duplicated
      
      #Variables to match
      mod_match <- mod_res[,c("experiment", "ensemble", "variable")]
      
      #Find duplicates
      vr_ind <- which(duplicated(mod_match))
      
      #If found dupcliates, only get latest version
      if (length(vr_ind) > 0) {
        
        rm_ind <- vector()
        
        for (v in 1:length(vr_ind)) {
          
          #Finds rows that are duplicates
          common_ind <- which(apply(mod_match, MARGIN=1, function(x) all(x == mod_match[vr_ind[v],])))
          
          #Check which one of these is the newest version
          #(remove "v" from start, convert to numeric and find biggest)
          max_ind <- which.max(as.numeric(substr(final_results$version[common_ind], 2, 9)))
          
          
          rm_ind <- append(rm_ind, common_ind[-max_ind])
          
        }
        
        #Remove all old versions
        mod_res <- mod_res[-rm_ind,]
        
      }
      
      
      
      ### Copy files ###
      
      #Not saving these separately for each experiment at this stage,
      #Should probably be changed later
      
      #If saving to same directory
      if (combine) {
        #Create output directory
        target_dir <- paste(outdir, "/../Land_masks/", mod_res$model[1], sep="/")
        
        #Else
      } else {
        #Create output directory
        target_dir <- paste(outdir,  "/../Land_masks/", mod_res$model[1], sep="/")
      }
      
      dir.create(target_dir, recursive=TRUE, showWarnings = FALSE)
      
      
      #Loop through files
      for (f in 1:nrow(mod_res)) {
        
        #Basic sanity checks:
        #Check that path contains correct model, variable and experiment  
        
        #Extract data for this iteration
        entry <- mod_res[f,]
        
        
        #Model check
        if (!grepl(entry$model, entry$path)) {
          stop("File path doesn't contain correct model")
        } 
        #Variable check
        if (!grepl(entry$variable, entry$path)) {
          stop("File path doesn't contain correct variable")
        }
        
        
        #Experiment check
        if (!grepl(entry$experiment, entry$path)) {
          stop("File path doesn't contain correct experiment")
        }
        
        
        
        #Create symbolic link to file (need to add pattern because some models 
        #will otherwise return surplus variables. Also adding "_" to avoid similar 
        #var names to be returned)
        
        file.symlink(from=list.files(entry$path, full.names=TRUE, 
                                     pattern=paste0(entry$variable, "_")), 
                     to=target_dir)
        
        
        
        
        
      }
      
      
      
      
    } #if found files 
    
  } #models
  
}  
