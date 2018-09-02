#! /home/ekalina/miniconda2/bin/python

import os, sys, ConfigParser

# Read the configuration file.
config = ConfigParser.SafeConfigParser(os.environ)
config.read(['../parm/bwrf.conf'])

# Add the scripts directory to the system path
# and import the modules.
SCRIPTSbwrf = config.get("DEFAULT","SCRIPTSbwrf")
sys.path.append(SCRIPTSbwrf)
import launcher, inputs, wps, forecast, post

# Run the bWRF system.
launcher.launch(config)
inputs.init_input(config)
inputs.run_input(config)
wps.init_wps(config)
wps.run_geogrid(config)
wps.run_ungrib(config)
wps.run_metgrid(config)
forecast.init_wrf(config)
forecast.run_real(config)
forecast.run_forecast(config)
post.init_post(config)
post.run_post(config)
