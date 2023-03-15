Searching and post-processing CMIP data
=======================================

Searches NCI database for CMIP data to match models and ensembles across experiments, using most recent file versions. Then post-processes data (unit conversions, land masking, swaps hemispheres)


To process **CMIP6**, run:
`Process_CMIP6_outputs.sh`


To process **CMIP5**, run:
`Process_CMIP5_outputs.sh`

In either code, set options under "1. SET PATHS" and "2. SET OPTIONS". Then qsub script and you're off!


To see instructions for setting up the Clef search, run `clef cmip6 --help` for CMIP6 and `clef cmip5 --help` for CMIP5. 


Requires membership of NCI projects hh5 (for CMS python installation), and oi10 and fs38 (for CMIP data)

