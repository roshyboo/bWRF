#! python

import os, sys

# What is the scripts directory?
SCRIPTSbwrf=os.path.join(HOMEbwrf+'/scripts')

# Where is the parm directory?
PARMbwrf=os.path.join(HOMEbwrf+'/parm')

# Add the scripts directory to the system path.
sys.path.append(SCRIPTSbwrf)

# Import various packages from the scripts directory.
import ConfigParser
cp = ConfigParser.RawConfigParser()

import launch, inputter, wps, forecast, post

# Designate and create some directories.
CYCLEpath=launch.make_dirs(WORKbwrf+'/'+exp+'/'+YMDH+'/','-p')
INPUTpath=launch.make_dirs(CYCLEpath+'bwrfdata')
WPSpath=launch.make_dirs(CYCLEpath+'wps/')
UNGRIBpath=launch.make_dirs(WPSpath+'ungrib/')
METGRIBpath=launch.make_dirs(WPSpath+'metgrid/')
FORECASTpath=launch.make_dirs(CYCLEpath+'fcst/')
POSTpath=launch.make_dirs(CYCLEpath+'post/')

# Cat the confs.
CONFpath = CYCLEpath+exp+'/'+YMDH+".conf"
os.system("cat "+PARMbwrf+"/*.conf >& "+CONFpath)

# Read the cycle conf.
cp.read(CONFpath)
gfs_dataset=cp.get(inputter,gfs)
gfs_grib_item=cp.get(inputter,gfs_gribA)

inputter.fetch_item(gfs_url,gfs_grib_item,INPUTpath)
###wps.geogrid()
wps.ungrib.linkgrib(UNGRIBpath)
wps.ungrib.run_ungrib()
wps.metgrid.run_metgrid()
forecast.run_forecast()
post.run_post(FORECASTpath,POSTpath)
