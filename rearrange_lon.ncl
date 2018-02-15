load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("reorder_lon_dimTLLL")
function reorder_lon_dimTLLL(in_data[*][*][*][*]:numeric)
  local nhalf,nhalf_minus,dim_lengths, out_data,nsize
  begin

  dim_lengths = dimsizes(in_data)
  nsize = dim_lengths(3)
  nhalf = dim_lengths(3)/2
  nhalf_minus = nhalf - 1
  out_data = in_data
  out_data(:,:,:,0:nhalf_minus) = in_data(:,:,:,nhalf:nsize-1)
  out_data(:,:,:,nhalf:nsize-1) = in_data(:,:,:,0:nhalf-1)
  return(out_data)
end

undef("reorder_lon_dimTLL")
function reorder_lon_dimTLL(in_data[*][*][*]:numeric)
  local nhalf,nhalf_minus,dim_lengths, out_data,nsize
  begin

  dim_lengths = dimsizes(in_data)
  nsize = dim_lengths(2)
  nhalf = dim_lengths(2)/2
  nhalf_minus = nhalf - 1
  out_data = in_data
  out_data(:,:,0:nhalf_minus) = in_data(:,:,nhalf:nsize-1)
  out_data(:,:,nhalf:nsize-1) = in_data(:,:,0:nhalf-1)
  return(out_data)
end

undef("reorder_lon_dimLL")
function reorder_lon_dimLL(in_data[*][*]:numeric)
  local nhalf,nhalf_minus,dim_lengths, out_data,nsize
  begin

  dim_lengths = dimsizes(in_data)
  nsize = dim_lengths(1)
  nhalf = dim_lengths(1)/2
  nhalf_minus = nhalf - 1
  out_data = in_data
  out_data(:,0:nhalf_minus) = in_data(:,nhalf:nsize-1)
  out_data(:,nhalf:nsize-1) = in_data(:,0:nhalf-1)

  return(out_data)
end

undef("reorder_lon_dimTL")
function reorder_lon_dimTL(in_data[*][*]:numeric)
  local nhalf,nhalf_minus,dim_lengths, out_data,nsize
  begin

  dim_lengths = dimsizes(in_data)
  nsize = dim_lengths(1)
  nhalf = dim_lengths(1)/2
  nhalf_minus = nhalf - 1
  out_data = in_data
  out_data(:,0:nhalf_minus) = in_data(:,nhalf:nsize-1)
  out_data(:,nhalf:nsize-1) = in_data(:,0:nhalf-1)

  return(out_data)
end

undef("reorder_lon_dimLT")
function reorder_lon_dimLT(in_data[*][*]:numeric)
  local nhalf,nhalf_minus,dim_lengths, out_data,nsize
  begin

  dim_lengths = dimsizes(in_data)
  nsize = dim_lengths(1)
  nhalf = dim_lengths(1)/2
  nhalf_minus = nhalf - 1
  out_data = in_data
  out_data(0:nhalf_minus,:) = in_data(nhalf:nsize-1,:)
  out_data(nhalf:nsize-1,:) = in_data(0:nhalf-1,:)

  return(out_data)
end

undef("reorder_lon_dim")
function reorder_lon_dim(in_data[*]:numeric)
  local nhalf,nhalf_minus,dim_lengths, out_data,nsize
  begin

  dim_lengths = dimsizes(in_data)
  nsize = dim_lengths(0)
  nhalf = dim_lengths(0)/2
  nhalf_minus = nhalf - 1
  out_data = in_data
  out_data(0:nhalf_minus) = in_data(nhalf:nsize-1)
  out_data(nhalf:nsize-1) = in_data(0:nhalf-1)
  return(out_data)
end


undef("change_lon_range")
function change_lon_range(in_data[*]:numeric)
local nhalf,nhalf_minus,dim_lengths, out_data,nsize
begin

dim_lengths = dimsizes(in_data)
nsize = dim_lengths(0)
nhalf = dim_lengths(0)/2
nhalf_minus = nhalf - 1
out_data = in_data
out_data(0:nhalf_minus) = in_data(nhalf:nsize-1)
out_data(nhalf:nsize-1) = in_data(0:nhalf-1)
return(out_data)
end
