"""

GRID repair

Takes a netcdf file name as argument and fixes the lon grid


"""
#==============================================================================
__title__ = "Netcdf CMIP5 grid repair"
__author__ = "Arden Burrell"
__version__ = "1.0 (04.04.2018)"
__email__ = "arden.burrell@gmail.com"

#==============================================================================
# Import packages
import numpy as np 
import pandas as pd
import pdb
import sys
import os
import subprocess as subp
import argparse
import warnings as warn
#==============================================================================

def main(args):
	# ========== Load in the lookup table of the CMIP5 grids ==========
	# SOURCE: https://portal.enes.org/data/enes-model-data/cmip5/resolution
	grids   = pd.read_csv("./CMIP5model_lookup.csv", header=0)

	# ========== Set up the Key Infomation ==========
	fname   = args.fname
	# CHeck if the file ends with .nc
	if fname.endswith(".nc"):
		fname = fname[:-3]

	# ========== Create a for loop to loop over every model ==========
	fcleanup = repair_netcdf(fname, grids)

	# ========== cleanup the files ==========
	for file in fcleanup:
		os.remove(file)


#==============================================================================
def repair_netcdf(fname, grids):
	"""
	Repairs the grid of a given netcdf file
	"""

	# ========== Set the path and the file name ==========
	# fname = "%s_%s_%s_r1i1p1_%s_1950_2050_%s_regrid.nc" %(var, model, sen, units, sen)
	fout  = "%s_setgrid" % (fname)
	
	# ========== Create a list of files to cleanup ==========
	cleanup = []

	# ========== Check if the file exists ==========
	if not os.path.isfile(fname+".nc"):
		# check if the file exists with a different name
		raise IOError("WARNING: The file %s cannot be found"% fname)

	# ========== Create a new grid ==========
	# Save the current grid
	subp.call("cdo griddes %s.nc > %sGriddes" % (fname, fname), shell=True)
	# add the griddes to the cleanup 
	cleanup.append("%sGriddes" % fname)

	# open the current grid
	gfile    = open("%sGriddes" % fname, "r") 
	# Split the lines of the grid file
	ginfo    =  gfile.read().splitlines()
	# Check and see if the start is known
	if (
		any([n.startswith("xfirst") for n in ginfo])
		) and (
		any([n.startswith("xinc") for n in ginfo])
		):
		addxdet = False
		# Set the lines to be removed
		badel    = ["xvals", "yvals", "     ", "xbounds", "ybounds"]
	else:
		addxdet = True
		# Set the lines to be removed
		badel    = ["xvals", "yvals", "     ", "xbounds", "ybounds", "xfirst", "xinc"]

	# Create list to hold the new grid details
	new_grid = []

	for ginf in ginfo:
		test = []
		for be in badel:
			if ginf.startswith(be):
				test.append(False)
			elif ginf == "#":
				test.append(False)
			else:
				test.append(True)
		
		if all(test):
			new_grid.append(ginf)
	# Add the additional x variables
	if addxdet:
		# work out the model from the fname
		model = fname.split("/")[-2]
		new_grid.append('xfirst    = -180')
		new_grid.append('xinc      = %s' %  str(
			float(grids[grids["Model"]==model]["Longitude"]) ))
	

	# Check the y values, if they are missing use the ones in the original grid file
	if not (any([n.startswith("yfirst") for n in ginfo])):
		# print ("Seting the y bounds")
		vals = []
		for glov in range(0,len(ginfo)):
			if  ginfo[glov].startswith("yvals"):
				vals.append(glov)
			elif ginfo[glov].startswith("ybounds"):
				vals.append(glov)
		if len (vals) == 2:
			for yv in ginfo[vals[0]:vals[1]]:
				new_grid.append(yv)

		else:
			print("\n")
			raise IndexError("Bounding is incorrect")

	# Save the grid out
	newgrid = save_grid(fname, new_grid)
	cleanup.append(newgrid)

	# ========== Set the new grid file ==========
	# Save the current grid
	subp.call("cdo setgrid,%sGridFix %s.nc %s.nc" % (fname, fname, fout), shell=True)
	
	if not os.path.isfile("%s.nc" % fout):
		raise IOError("The output file was not created, going interactive")
	
	# ========== return the files to be removed ==========
	cleanup.append("%s.nc" % fname)
	return cleanup

#==============================================================================

def save_grid(fname, grid):
	"""Takes a list of elements and save them too a grid"""
	with open(("%sGridFix" % fname), 'w') as file_handler:
	    for item in grid:
	        file_handler.write("{}\n".format(item))
	        pass
    # Return the name of the file
	return ("%sGridFix" % fname)

#==============================================================================

if __name__ == '__main__':
	description='Arguments for grid repair'
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument(
		'fname', type=str, 
		help='The fname of the netcdf to have its grid repaired')
	args = parser.parse_args() 
	main(args)



