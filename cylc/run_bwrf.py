#! /home/ekalina/miniconda2/bin/python

import os, sys

# Get arguments passed to script.
#YMDH=sys.argv[1]
#expt=sys.argv[2] 
#HOMEbwrf=sys.argv[3]
#WORKbwrf=sys.argv[4]

# Make a unique work directory for this cycle.
#WORKbwrf=WORKbwrf+'/'+expt+'/'+YMDH+'/'

# What is the scripts directory?
#SCRIPTSbwrf=os.path.join(HOMEbwrf+'/scripts')

# Where is the parm directory?
#PARMbwrf=os.path.join(HOMEbwrf+'/parm/')

# Import various packages from the scripts directory.
import ConfigParser
config = ConfigParser.SafeConfigParser(os.environ)
config.read(['../parm/bwrf.conf'])

# Add the scripts directory to the system path.
SCRIPTSbwrf = config.get("DEFAULT","SCRIPTSbwrf")
sys.path.append(SCRIPTSbwrf)

import inputs, wps, forecast, post

# Designate and create some directories.
#INPUTpath=WORKbwrf+'bwrfdata/'
#WPSpath=WORKbwrf+'wps/'
#UNGRIBpath=WPSpath+'ungrib/'
#METGRIBpath=WPSpath+'metgrid/'
#FORECASTpath=WORKbwrf+'fcst/'
#POSTpath=WORKbwrf+'post/'

#os.system("mkdir -p "+WORKbwrf)

#os.system("mkdir "+INPUTpath)
#os.system("mkdir "+WPSpath)
#os.system("mkdir "+UNGRIBpath)
#os.system("mkdir "+METGRIBpath)
#os.system("mkdir "+FORECASTpath)
#os.system("mkdir "+POSTpath)

# Cat the confs.
#os.system("cat "+PARMbwrf+"/*.conf > "+WORKbwrf+YMDH+".conf")

# Read the cycle conf.
#cp.read(WORKbwrf+YMDH+".conf")
gfs_dataset=config.get("inputs","gfs")
gfs_grib_item=config.get("inputs","gfs_gribA")

# Expand the gfs grib item into a list of files.
run_hours = config.getfloat("wrf","run_hours")
interval_hours = config.getfloat("wrf","interval_seconds")/3600

ihr = 0
while ihr <= run_hours:
  fahr = int(ihr)
  config.set("inputs","fahr",'{:03d}'.format(fahr))
  gfs_grib_item = config.get("inputs","gfs_gribA")
#  inputs.fetch_item(gfs_dataset,gfs_grib_item,config)
  ihr = ihr + interval_hours

#wps.init_wps(config)
#wps.run_geogrid(config)
#wps.run_ungrib(config)
#wps.run_metgrid(config)
forecast.init_wrf(config)
#forecast.run_real(config)
forecast.run_forecast(config)
#post.run_post(FORECASTpath,POSTpath)
