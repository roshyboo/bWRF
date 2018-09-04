#! /home/ekalina/miniconda2/bin/python

import os, sys, ConfigParser, time

start = time.time()

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

input_start=time.time()
inputs.init_input(config)
inputs.run_input(config)
print("Total time in input: ", (time.time()-input_start)/60.0, " minutes")

wps_start=time.time()
wps.init_wps(config)
wps.run_geogrid(config)
wps.run_ungrib(config)
wps.run_metgrid(config)
print("Total time in wps: ", (time.time()-wps_start)/60.0, " minutes")

fcst_start=time.time()
forecast.init_wrf(config)
forecast.run_real(config)
forecast.run_forecast(config)
print("Total time in fcst: ", (time.time()-fcst_start)/60.0, " minutes")

post_start=time.time()
post.init_post(config)
post.run_post(config)
print("Total time in post: ", (time.time()-post_start)/60.0, " minutes")

print("Total time to run bWRF: ", (time.time()-start)/60.0, " minutes")
