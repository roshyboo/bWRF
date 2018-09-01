#! /home/ekalina/miniconda2/bin/python

import os, sys, datetime, ConfigParser, f90nml

# Read the configuration file.
config = ConfigParser.SafeConfigParser(os.environ)
config.read(['../parm/bwrf.conf'])

# Add the scripts directory to the system path
# and import the modules.
SCRIPTSbwrf = config.get("DEFAULT","SCRIPTSbwrf")
sys.path.append(SCRIPTSbwrf)
import inputs, wps, forecast, post

# Is this an operational run?
# If so, figure out the cycle to run.
aYMDH = config.get("DEFAULT","aYMDH")
aHH = config.get("DEFAULT","aHH")

if aYMDH == "AUTO":
  currentYMD = datetime.datetime.now()
  aYMD = currentYMD.strftime("%Y%m%d")
  aYMDH = aYMD + aHH
  print("Will run the "+aYMDH+" cycle")
  config.set("DEFAULT","aYMDH",aYMDH)

# Create the work directory.
WORKbwrf=config.get("DEFAULT","WORKbwrf")

if not os.path.exists(WORKbwrf):
  os.makedirs(WORKbwrf)

# Update the namelists for this configuration.
wps_nml_file = config.get("wps","namelist")
wrf_nml_file = config.get("wrf","namelist")
WPShome = config.get("DEFAULT","WPShome")
WRFhome = config.get("DEFAULT","WRFhome")
wps_nml_file = config.get("wps","namelist")
wrf_nml_file = config.get("wrf","namelist")
run_hours = config.getint("wrf","run_hours")

wps_nml = f90nml.read(WPShome+"/"+wps_nml_file)
wrf_nml = f90nml.read(WRFhome+"/"+wrf_nml_file)

aYMDH_dt = datetime.datetime.strptime(aYMDH,'%Y%m%d%H')
wps_nml['share']['start_date'] = aYMDH_dt.strftime('%Y-%m-%d_%H:00:00')
wps_nml['share']['end_date'] = (aYMDH_dt + datetime.timedelta(hours=run_hours)).strftime('%Y-%m-%d_%H:00:00')

wrf_nml['time_control']['start_year'] = '{:04d}'.format(int(aYMDH[0:4]))
wrf_nml['time_control']['start_month'] = '{:02d}'.format(int(aYMDH[4:6]))
wrf_nml['time_control']['start_day'] = '{:02d}'.format(int(aYMDH[6:8]))
wrf_nml['time_control']['start_hour'] = '{:02d}'.format(int(aYMDH[8:10]))

wps_nml.write(WPShome+"/"+wps_nml_file, force=True)
wrf_nml.write(WRFhome+"/"+wrf_nml_file, force=True)

# Run the bWRF system.
#inputs.init_input(config)
#inputs.run_input(config)
#wps.init_wps(config)
#wps.run_geogrid(config)
#wps.run_ungrib(config)
#wps.run_metgrid(config)
#forecast.init_wrf(config)
#forecast.run_real(config)
#forecast.run_forecast(config)
#post.init_post(config)
#post.run_post(config)
