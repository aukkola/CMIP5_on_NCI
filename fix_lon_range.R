fix_lon_range <- function(file){
  
  library(ncdf4)
  
  #Open file handle
  nc <- nc_open(file, write=TRUE,readunlim=FALSE)
  
  #Get longitude variables
  lon      <- ncvar_get(nc, "lon")
  
  #If already fixed, skip
  if(min(lon) < -170){
    return()
  }
    
  #Correct lon
  lon <- lon - 180
  
  #Get lon bounds (try both varnames, varies by file)
  lon_bnds <- tryCatch(ncvar_get(nc, "lon_bounds") - 180, 
                       error = function(e) "try-error")

    varname  <- "lon_bounds"
  
  if(lon_bnds == "try-error"){
    lon_bnds <- ncvar_get(nc, "lon_bnds") - 180
    varname  <- "lon_bnds"
  }
  
  #Save new variables to file
  ncvar_put(nc, "lon", lon)
  ncvar_put(nc, varname, lon_bnds)
  
  #Close file
  nc_close(nc)

}



