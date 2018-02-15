load "$NCARG_LIB/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_LIB/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_LIB/ncarg/nclscripts/csm/contributed.ncl"
load "rearrange_lon.ncl"

  setvalues NhlGetWorkspaceObjectId()
    "wsMaximumSize" : 1000000000
  end setvalues


begin

possible_lon_names = (/"lon","longitude","LON","LONGITUDE","x","X"/)


file_input = getenv("INPUTFILE")

file_output = getenv("OUTPUTFILE")

nc_in = addfile(file_input,"r")

system("rm "+file_output)

nc_out = addfile(file_output,"c")


var_names = getfilevarnames(nc_in)
nvars = dimsizes(var_names)

do ivar=0,nvars-1

   tmp_name = var_names(ivar)
   tmp = nc_in->$tmp_name$

   tmp_att = getvaratts(tmp)

   tmp_dims = getvardims(tmp)

   ndims = dimsizes(tmp_dims)

   reorder = False

   do idim=0,dimsizes(possible_lon_names)-1
      if (isdim(tmp,possible_lon_names(idim))) then
         actual_name = possible_lon_names(idim)
         reorder = True
      end if
   end do

  printVarSummary(tmp)

   if (reorder .eq. True) then
   ;find the dim I need to reorder
      do idim=0,ndims-1
         if (.not.ismissing(tmp_dims(idim))) then
            if (tmp_dims(idim) .eq. actual_name) then
              actual_idim = idim
            end if
         end if
      end do

      print(actual_idim)
      print(dimsizes(tmp_dims))

      if ((actual_idim .eq. 0) .and. (dimsizes(tmp_dims) .eq. 1)) then
         tmp_new = reorder_lon_dim(tmp)
      end if

      if ((actual_idim .eq. 0) .and. (dimsizes(tmp_dims) .eq. 2)) then
         tmp_new = reorder_lon_dimLT(tmp)
      end if

      if ((actual_idim .eq. 1) .and. (dimsizes(tmp_dims) .eq. 2)) then
         tmp_new = reorder_lon_dimLL(tmp)
      end if

      if ((actual_idim .eq. 2) .and. (dimsizes(tmp_dims) .eq. 3)) then
         tmp_new = reorder_lon_dimTLL(tmp)
      end if

      if ((actual_idim .eq. 3) .and. (dimsizes(tmp_dims) .eq. 4)) then
         tmp_new = reorder_lon_dimTLLL(tmp)
      end if

      copy_VarAtts(tmp,tmp_new)

      nc_out->$tmp_name$=tmp_new

      ;Change longitude range back so matches lon_bounds (goes 0-360 which is not correct but at least consistent)
      ;Not yet sure how to adjust lon_bounds
	  lon_temp = change_lon_range(nc_out&$actual_name$)
	  nc_out&$actual_name$=lon_temp


      delete(tmp_new)


      else
      	nc_out->$tmp_name$=nc_in->$tmp_name$
      end if


   delete(tmp)
   delete(tmp_att)
   delete(tmp_dims)

end do


delete(nc_out)
delete(nc_in)


end


