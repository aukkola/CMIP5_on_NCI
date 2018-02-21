### Find and sym-link copy CMIP5 runs to match criteria ###

library(DBI)
library(RSQLite)


#clear R environment
rm(list=ls(all=TRUE))


####################
### Set criteria ###
####################

### 1. SET PATH ###

#Directory where to copy data
outdir <- "/g/data1/w35/amu561/CMIP5_fluxnet/CMIP5_Data"


### 2. SELECT VARIABLES AND EXPERIMENTS ###

#Variables and experiments
variables  <- c("pr", "tas")
experiment <- c("historical", "rcp85")
mip        <- c("Amon")


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
combine  <- TRUE
dir_name <- "historical_rcp8.5" 


### 5. DECIDE IF WANT LAND MASKS ###

#Retrieves land masks for selected models 
get_land_masks <- TRUE



#------------------------------------------------------------------------------

#Connect to database (this points to the latest database version)
db  <- "/g/data1/ua6/unofficial-ESG-replica/tmp/tree/cmip5_raijin_latest.db"
con <- dbConnect(SQLite(), dbname=db)
on.exit(dbDisconnect(con))




######################
### Query dataset ####
######################

#Create query string
#Need to use paste when using R variables in query

#Create strings for search
exp_str <- paste(paste("'", experiment, "'", sep=""), collapse=",")
mip_str <- paste(paste("'", mip, "'", sep=""), collapse=",")
var_str <- paste(paste("'", variables, "'", sep=""), collapse=",")
ens_str <- paste(paste("'", ensemble, "'", sep=""), collapse=",")



#Create complete query string
if (!is.na(ensemble[1])) {
  
  query <- paste("SELECT * FROM instances WHERE experiment in (", exp_str, 
                 ") and mip in (", mip_str, ") and variable in (", var_str, ")",
                 " and ensemble in (", ens_str, ")", sep="")
  
} else {
  
  query <- paste("SELECT * FROM instances WHERE experiment in (", exp_str, 
                 ") and mip in (", mip_str, ") and variable in (", var_str, ")",
                 sep="")
}


#Perform query
results <- dbGetQuery(con, query)





##########################
### Find common models ###
##########################

#Separate results by experiment
res_experiment <- lapply(experiment, function(x) results[which(results$experiment==x),])


#Find models for each experiment and variable
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


selected_ens <- vector()
for (k in 1:length(common_models)) {
  
  #Find all ensemble members for each experiment and variable
  avail_ens <- lapply(1:length(ensemble), function(x) mapply(function(y,z) y[which(z==common_models[k])], 
                                                             y=ensemble[[x]], z=models[[x]], SIMPLIFY=FALSE))
  
  
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
  
  #Select first available ensemble
  if (length(common_ens) > 0) {
    #Sorting this way so r1i1p1 comes before r10i1p1 etc.
    selected_ens[k] <- common_ens[order(nchar(common_ens), common_ens)][1]
  }
}


# Check if didn't find any common ensemble members for a model
# Remove model in this case
if (any(is.na(selected_ens))) {
  common_models <- common_models[-which(is.na(selected_ens))]
  selected_ens <- selected_ens[-which(is.na(selected_ens))]
}


if (exists("mods")) {
  if (!is.na(mods[1])) {
    
    #find indices for required models
    req_models <- sapply(mods, function(x) which(common_models==x))
    
    common_models <- common_models[req_models]
    selected_ens  <- selected_ens[req_models]
    
  }
}


#Find instance numbers for selected models and ensembles

instances <- lapply(res_experiment, function(x) mapply(function(y,z) x$instance_id[which(x$model==y & x$ensemble==z)], 
                                                       y=common_models, z=selected_ens))

instances <- unlist(instances)


### Find paths for selected models and ensembles ###

#Create string
inst_str <- paste("(", paste(instances, collapse=","), ")", sep="")
query_final <- paste("SELECT * FROM versions WHERE instance_id in", inst_str)

paths <- dbGetQuery(con, query_final)





############################
### Select file versions ###
############################


versions <- vector()

for (k in 1:length(instances)) {
  
  #Find indices for instance
  ind <- which(paths$instance_id==instances[k])
  
  #Extract all results
  all_versions <- paths[ind,]
  
  #If only one version, use that
  if (nrow(all_versions)==1) {
    
    versions[k] <- all_versions$version_id
    
  #If multiple versions, find newest
  #NB. sometimes latest version does not contain all the files
  #e.g. HadGEM2-ES day/pr. Not implemented a fix for this yet
    
  } else {
  
    #Find latest  
    latest <- which(all_versions$is_latest==1)
    
    
    #Found one latest, use that
    if (length(latest)==1) {
      
      versions[k] <- all_versions$version_id[latest]
      
    #Found more than one latest, compare versions
    } else if (length(latest)>1) {
      
      #If all versions the same, use first result
      if (all(all_versions$version[latest] == all_versions$version[latest[1]])) {  
        versions[k] <- all_versions$version_id[latest[1]]
        
      } else {
        
        #Sort and pick first (=newest)
        
        stop("This bit of code needs FIXING")
        sorted <- sort(all_versions$version[latest], decreasing=TRUE)
        
        versions[k] <- all_versions$version[sorted[1]]
      }
      
      
    #Didn't find a latest, sort versions and pick latest
    } else {  
      
      #Sort decreasing
      sorted <- sort(all_versions$version, decreasing=TRUE)
      
      #Pick first instance and find corresponding version no.
      versions[k] <- all_versions$version_id[which(all_versions$version == sorted[1])[1]]
    
    }
  }
}





### Collate final search results ###

#Find versions to extract
ind <- sapply(versions, function(x) which(paths$version_id==x))

#Extract final paths
final_paths <- paths[ind,]
  
#Find instance numbers
final_instances <- paths$instance_id[ind]

#Extract final models
ind1 <- sapply(final_instances, function(x) which(results$instance_id==x))
final_models <- results[ind1,]



###########################################
### Create database of selected outputs ###
###########################################



for (k in 1:nrow(final_models)) {
  
  #If saving to same directory
  if (combine) {
    #Create output directory
    target_dir <- paste(outdir, dir_name, final_models$mip[k],
                        final_models$variable[k], final_models$model[k],
                        final_models$ensemble[k], sep="/")

  #Else
  } else {
    #Create output directory
    target_dir <- paste(outdir, final_models$experiment[k], final_models$mip[k],
                        final_models$variable[k], final_models$model[k],
                        final_models$ensemble[k], sep="/")
  }
  
  dir.create(target_dir, recursive=TRUE, showWarnings = FALSE)
  
  
  
  
  #Basic sanity checks:
  #Check that path contains correct model, variable and experiment  
  
  #Need an exception for the GISS and FIO-ESM models as file path structure is different
  
  #GISS
  if (any(final_models$model[k]==c("GISS-E2-R", "GISS-E2-H", "GISS-E2-H-CC", "GISS-E2-R-CC"))) {
    #Get model version
    ext <- strsplit(final_models$model[k], "GISS-")[[1]][2]
    #Check
    if (!grepl("GISS", final_paths$path[k]) & !grepl(ext, final_paths$path[k])) {
      stop("File path doesn't contain correct model")
    }
    
  #FIO-ESM (slightly different spelling in paths)  
  } else if (final_models$model[k] == "FIO-ESM") {
    if (!grepl("FIO_ESM", final_paths$path[k])) {
      stop("File path doesn't contain correct model")
    }  
  
  #Other models
  } else {
    #Model check
    if (!grepl(final_models$model[k], final_paths$path[k])) {
    stop("File path doesn't contain correct model")
    } 
    #Variable check
    if (!grepl(final_models$variable[k], final_paths$path[k])) {
      stop("File path doesn't contain correct variable")
    }
  }  
  
  #Experiment check
  if (!grepl(final_models$experiment[k], final_paths$path[k])) {
    stop("File path doesn't contain correct experiment")
  }



  
  #Create symbolic link to file (need to add pattern because some models 
  #will otherwise return surplus variables. Also adding "_" to avoid similar 
  #var names to be returned)

  file.symlink(from=list.files(final_paths$path[k], full.names=TRUE, 
                               pattern=paste0(final_models$variable[k], "_")), 
               to=target_dir)
  
}




#########################################
### Retrieve corresponding land masks ###
#########################################


if (get_land_masks) {
  
  #Mask variable name
  mask_var  <- "sftlf"
  
  #Create string for selected models
  mod_str <- paste(paste("'", unique(final_models$model), "'", sep=""), collapse=",")
  
  
  #Create query string (only use first experiment if combining experiments)
  if (combine) {
    mask_query <- paste("SELECT * FROM instances WHERE experiment in ('", experiment[1], 
                        "') and variable in ('", mask_var, "') and model in (",
                        mod_str, ") and ensemble in ('r0i0p0')", sep="")
    
  } else {
    mask_query <- paste("SELECT * FROM instances WHERE experiment in (", exp_str, 
                        ") and variable in ('", mask_var, "') and model in (",
                        mod_str, ") and ensemble in ('r0i0p0')", sep="")
  }
  
  
  #Perform query
  mask_results <- dbGetQuery(con, mask_query)
  
  
  
  ### Find paths for masks ###
  
  instances <- mask_results$instance_id
  
  #Create string
  inst_str    <- paste("(", paste(instances, collapse=","), ")", sep="")
  query_final <- paste("SELECT * FROM versions WHERE instance_id in", inst_str)
  
  mask_paths <- dbGetQuery(con, query_final)
  
  
  
  ############################
  ### Select file versions ###
  ############################
  
  
  versions <- vector()
  
  for (k in 1:length(instances)) {
    
    #Find indices for instance
    ind <- which(mask_paths$instance_id==instances[k])
    
    #Extract all results
    all_versions <- mask_paths[ind,]
    
    #If only one version, use that
    if (nrow(all_versions)==1) {
      
      versions[k] <- all_versions$version_id
      
      #If multiple versions, find newest
    } else {
      
      #Find latest  
      latest <- which(all_versions$is_latest==1)
      
      
      #Found one latest, use that
      if (length(latest)==1) {
        
        versions[k] <- all_versions$version_id[latest]
        
        #Found more than one latest, compare versions
      } else if (length(latest)>1) {
        
        #If all versions the same, use first result
        if(all(all_versions$version[latest] == all_versions$version[latest[1]])){  
          versions[k] <- all_versions$version_id[latest[1]]
          
        } else {
          
          #Sort and pick first (=newest)
          
          stop("This bit of code needs FIXING")
          sorted <- sort(all_versions$version[latest], decreasing=TRUE)
          
          versions[k] <- all_versions$version[sorted[1]]
        }
        
        
        #Didn't find a latest, sort versions and pick latest
      } else {  
        
        #Sort decreasing
        sorted <- sort(all_versions$version, decreasing=TRUE)
        
        #Pick first instance and find corresponding version no.
        versions[k] <- all_versions$version_id[which(all_versions$version == sorted[1])[1]]
        
      }
    }
  }
  
  
  ### Collate final search results ###
  
  #Find versions to extract
  ind <- sapply(versions, function(x) which(mask_paths$version_id==x))
  
  #Extract final paths
  final_paths <- mask_paths[ind,]
  
  #Find instance numbers
  final_instances <- mask_paths$instance_id[ind]
  
  #Extract final models
  ind1 <- sapply(final_instances, function(x) which(mask_results$instance_id==x))
  final_models_mask <- mask_results[ind1,]
  
  
  #Check if found a mask for all models, return a warning if not
  common_mods <- intersect(final_models$model, final_models_mask$model)
  
  if (length(common_mods) != length(unique(final_models$model))) {
    
   not_found <- is.element(unique(final_models$model), common_mods)
   
    warning(paste("Could not find masks for models: ", 
                  paste(unique(final_models$model)[!not_found],
                  collapse=", ")))
  }
  
    
  
  ###########################################
  ### Create database of selected outputs ###
  ###########################################
  
  
  
  for (k in 1:nrow(final_models_mask)) {
    
    #If saving to same directory
    if (combine) {
      #Create output directory
      target_dir <- paste(outdir, "../Processed_masks", dir_name,
                          final_models_mask$model[k], sep="/")
      
    #Else
    } else {
      #Create output directory
      target_dir <- paste(outdir,  "../Processed_masks", final_models_mask$experiment[k], 
                          final_models_mask$model[k], sep="/")
    }
    
    dir.create(target_dir, recursive=TRUE, showWarnings = FALSE)
    
    
    
    #Basic sanity checks:
    #Check that path contains correct model, variable and experiment  
    
    #Need an exception for the GISS and FIO-ESM models as file path structure is different
    
    #GISS
    if (any(final_models_mask$model[k]==c("GISS-E2-R", "GISS-E2-H", "GISS-E2-H-CC", "GISS-E2-R-CC"))) {
      #Get model version
      ext <- strsplit(final_models_mask$model[k], "GISS-")[[1]][2]
      #Check
      if (!grepl(ext, final_paths$path[k])) {
        stop("File path doesn't contain correct model")
      }

      #FIO-ESM (slightly different spelling in paths)
    } else if (final_models_mask$model[k] == "FIO-ESM") {
      if (!grepl("FIO_ESM", final_paths$path[k])) {
        stop("File path doesn't contain correct model")
      }

      #Other models
    } else {
      #Model check
      if (!grepl(final_models_mask$model[k], final_paths$path[k])) {
        stop("File path doesn't contain correct model")
      } 
      #Variable check
      if (!grepl(final_models_mask$variable[k], final_paths$path[k])) {
        stop("File path doesn't contain correct variable")
      }
    }  
    
    #Experiment check
    if (!grepl(final_models_mask$experiment[k], final_paths$path[k])) {
      stop("File path doesn't contain correct experiment")
    }
    
    
    
    
    #Create symbolic link to file (need to add pattern because some models 
    #will otherwise return surplus variables. Also adding "_" to avoid similar 
    #var names to be returned)
    
    file.symlink(from=list.files(final_paths$path[k], full.names=TRUE, 
                                 pattern=paste0(final_models_mask$variable[k], "_")), 
                 to=target_dir)
    
  }
  
  
  
} #masks





