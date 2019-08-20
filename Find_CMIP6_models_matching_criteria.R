### Find and sym-link copy CMIP5 runs to match criteria ###


#clear R environment
rm(list=ls(all=TRUE))


####################
### Set criteria ###
####################

### 1. SET PATH ###

#Directory where to copy data
outdir <- "/g/data1/w35/amu561/CMIP6_drought/CMIP6_Data/Raw_CMIP6_data"


### 3. SELECT ENSEMBLES ###

#Select ensemble (set to NA if don't require a specific ensemble member,
#in this case code will select the first member common to all variables. If
#no common member found, model won't be processed)
ensemble <- NA #"r1i1p1"


### 4. DECIDE FILE STRUCTURE ###

#Should model files for all experiments be
#saved in same folder (specify name of folder)?
#If want to e.g. combine historical and RCP8.5 runs,
#use this option, else set to FALSE
combine  <- FALSE
dir_name <- "historical_rcp4.5" 


### 5. DECIDE IF WANT LAND MASKS ###

#Retrieves land masks for selected models 
get_land_masks <- TRUE



#------------------------------------------------------------------------------


##########################
### Get query dataset ####
##########################


#Get Clef search results
results <- as.vector(read.csv("/g/data1/w35/amu561/CMIP6_drought/CMIP6_Data/temp_res/cmip6_clef_search_results.csv", 
                    header=FALSE, colClasses="character"))



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
sorted_results$ensemble   <- sapply(all_res, function(x) x[11])

#Get variable
sorted_results$variable   <- sapply(all_res, function(x) x[13])

#Get grid type
sorted_results$grid       <- sapply(all_res, function(x) x[14])

#Get version
sorted_results$version    <- sapply(all_res, function(x) x[15])

#Get time resolution
sorted_results$time_resolution <- sapply(all_res, function(x) x[12])

#Get paths
sorted_results$path <- results[,1]



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





#############################################
### Should add a check for grid type here ###
#############################################







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
  ind <- which(final_results$model == common_models[k])
  
  #Find extra ensemble members
  ens_ind <- which(!(final_results[ind,]$ensemble %in% selected_ens[[k]]))
  
  if (length(ens_ind) >0) final_results <- final_results[-ind[ens_ind]]
  
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
    max_ind <- which.max(as.numeric(substr(final_results$version[common_ind], 2, 9)))

      
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



# 
# #########################################
# ### Retrieve corresponding land masks ###
# #########################################
# 
# 
# if (get_land_masks) {
#   
#   #Mask variable name
#   mask_var  <- "sftlf"
#   
#   #Create string for selected models
#   mod_str <- paste(paste("'", unique(final_models$model), "'", sep=""), collapse=",")
#   
#   
#   #Create query string (only use first experiment if combining experiments)
#   if (combine) {
#     mask_query <- paste("SELECT * FROM instances WHERE experiment in ('", experiment[1], 
#                         "') and variable in ('", mask_var, "') and model in (",
#                         mod_str, ") and ensemble in ('r0i0p0')", sep="")
#     
#   } else {
#     mask_query <- paste("SELECT * FROM instances WHERE experiment in (", exp_str, 
#                         ") and variable in ('", mask_var, "') and model in (",
#                         mod_str, ") and ensemble in ('r0i0p0')", sep="")
#   }
#   
#   
#   #Perform query
#   mask_results <- dbGetQuery(con, mask_query)
#   
#   
#   
#   ### Find paths for masks ###
#   
#   instances <- mask_results$instance_id
#   
#   #Create string
#   inst_str    <- paste("(", paste(instances, collapse=","), ")", sep="")
#   query_final <- paste("SELECT * FROM versions WHERE instance_id in", inst_str)
#   
#   mask_paths <- dbGetQuery(con, query_final)
#   
#   
#   
#   ############################
#   ### Select file versions ###
#   ############################
#   
#   
#   versions <- vector()
#   
#   for (k in 1:length(instances)) {
#     
#     #Find indices for instance
#     ind <- which(mask_paths$instance_id==instances[k])
#     
#     #Extract all results
#     all_versions <- mask_paths[ind,]
#     
#     #If only one version, use that
#     if (nrow(all_versions)==1) {
#       
#       versions[k] <- all_versions$version_id
#       
#       #If multiple versions, find newest
#     } else {
#       
#       #Find latest  
#       latest <- which(all_versions$is_latest==1)
#       
#       
#       #Found one latest, use that
#       if (length(latest)==1) {
#         
#         versions[k] <- all_versions$version_id[latest]
#         
#         #Found more than one latest, compare versions
#       } else if (length(latest)>1) {
#         
#         #If all versions the same, use first result
#         if(all(all_versions$version[latest] == all_versions$version[latest[1]])){  
#           versions[k] <- all_versions$version_id[latest[1]]
#           
#         } else {
#           
#           #Sort and pick first (=newest)
#           
#           stop("This bit of code needs FIXING")
#           sorted <- sort(all_versions$version[latest], decreasing=TRUE)
#           
#           versions[k] <- all_versions$version[sorted[1]]
#         }
#         
#         
#         #Didn't find a latest, sort versions and pick latest
#       } else {  
#         
#         #Sort decreasing
#         sorted <- sort(all_versions$version, decreasing=TRUE)
#         
#         #Pick first instance and find corresponding version no.
#         versions[k] <- all_versions$version_id[which(all_versions$version == sorted[1])[1]]
#         
#       }
#     }
#   }
#   
#   
#   ### Collate final search results ###
#   
#   #Find versions to extract
#   ind <- sapply(versions, function(x) which(mask_paths$version_id==x))
#   
#   #Extract final paths
#   final_paths <- mask_paths[ind,]
#   
#   #Find instance numbers
#   final_instances <- mask_paths$instance_id[ind]
#   
#   #Extract final models
#   ind1 <- sapply(final_instances, function(x) which(mask_results$instance_id==x))
#   final_models_mask <- mask_results[ind1,]
#   
#   
#   #Check if found a mask for all models, return a warning if not
#   common_mods <- intersect(final_models$model, final_models_mask$model)
#   
#   if (length(common_mods) != length(unique(final_models$model))) {
#     
#    not_found <- is.element(unique(final_models$model), common_mods)
#    
#     warning(paste("Could not find masks for models: ", 
#                   paste(unique(final_models$model)[!not_found],
#                   collapse=", ")))
#   }
#   
#     
#   
#   ###########################################
#   ### Create database of selected outputs ###
#   ###########################################
#   
#   
#   
#   for (k in 1:nrow(final_models_mask)) {
#     
#     #If saving to same directory
#     if (combine) {
#       #Create output directory
#       target_dir <- paste(outdir, "../Processed_masks", dir_name,
#                           final_models_mask$model[k], sep="/")
#       
#     #Else
#     } else {
#       #Create output directory
#       target_dir <- paste(outdir,  "../Processed_masks", final_models_mask$experiment[k], 
#                           final_models_mask$model[k], sep="/")
#     }
#     
#     dir.create(target_dir, recursive=TRUE, showWarnings = FALSE)
#     
#     
#     
#     #Basic sanity checks:
#     #Check that path contains correct model, variable and experiment  
#     
#     #Need an exception for the GISS and FIO-ESM models as file path structure is different
#     
#     #GISS
#     if (any(final_models_mask$model[k]==c("GISS-E2-R", "GISS-E2-H", "GISS-E2-H-CC", "GISS-E2-R-CC"))) {
#       #Get model version
#       ext <- strsplit(final_models_mask$model[k], "GISS-")[[1]][2]
#       #Check
#       if (!grepl(ext, final_paths$path[k])) {
#         stop("File path doesn't contain correct model")
#       }
# 
#       #FIO-ESM (slightly different spelling in paths)
#     } else if (final_models_mask$model[k] == "FIO-ESM") {
#       if (!grepl("FIO_ESM", final_paths$path[k])) {
#         stop("File path doesn't contain correct model")
#       }
# 
#       #Other models
#     } else {
#       #Model check
#       if (!grepl(final_models_mask$model[k], final_paths$path[k])) {
#         stop("File path doesn't contain correct model")
#       } 
#       #Variable check
#       if (!grepl(final_models_mask$variable[k], final_paths$path[k])) {
#         stop("File path doesn't contain correct variable")
#       }
#     }  
#     
#     #Experiment check
#     if (!grepl(final_models_mask$experiment[k], final_paths$path[k])) {
#       stop("File path doesn't contain correct experiment")
#     }
#     
#     
#     
#     
#     #Create symbolic link to file (need to add pattern because some models 
#     #will otherwise return surplus variables. Also adding "_" to avoid similar 
#     #var names to be returned)
#     
#     file.symlink(from=list.files(final_paths$path[k], full.names=TRUE, 
#                                  pattern=paste0(final_models_mask$variable[k], "_")), 
#                  to=target_dir)
#     
#   }
#   
#   
#   
# } #masks
# 
# 
# 
# 
# 
