import os, datetime, f90nml

def launch(conf):

# Is this an operational run?
# If so, figure out the cycle to run.
  aYMDH = conf.get("DEFAULT","aYMDH")
  aHH = conf.get("DEFAULT","aHH")

  if aYMDH == "AUTO":
    currentYMD = datetime.datetime.now()
    aYMD = currentYMD.strftime("%Y%m%d")
    aYMDH = aYMD + aHH
    print("Will run the "+aYMDH+" cycle")
    conf.set("DEFAULT","aYMDH",aYMDH)

# Create the work directory.
  WORKbwrf=conf.get("DEFAULT","WORKbwrf")

  if not os.path.exists(WORKbwrf):
    os.makedirs(WORKbwrf)

# Update the namelists for this configuration.
  wps_nml_file = conf.get("wps","namelist")
  wrf_nml_file = conf.get("wrf","namelist")
  WPShome = conf.get("DEFAULT","WPShome")
  WRFhome = conf.get("DEFAULT","WRFhome")
  wps_nml_file = conf.get("wps","namelist")
  wrf_nml_file = conf.get("wrf","namelist")
  run_hours = conf.getint("wrf","run_hours")

  wps_nml = f90nml.read(WPShome+"/"+wps_nml_file)
  wrf_nml = f90nml.read(WRFhome+"/"+wrf_nml_file)

  fYMDH_start = datetime.datetime.strptime(aYMDH,'%Y%m%d%H')
  fYMDH_end = fYMDH_start + datetime.timedelta(hours=run_hours)

  wps_nml['share']['start_date'] = fYMDH_start.strftime('%Y-%m-%d_%H:00:00')
  wps_nml['share']['end_date'] = fYMDH_end.strftime('%Y-%m-%d_%H:00:00')

  wrf_nml['time_control']['start_year'] = fYMDH_start.year
  wrf_nml['time_control']['start_month'] = fYMDH_start.month
  wrf_nml['time_control']['start_day'] = fYMDH_start.day
  wrf_nml['time_control']['start_hour'] = fYMDH_start.hour

  wrf_nml['time_control']['end_year'] = fYMDH_end.year
  wrf_nml['time_control']['end_month'] = fYMDH_end.month
  wrf_nml['time_control']['end_day'] = fYMDH_end.day
  wrf_nml['time_control']['end_hour'] = fYMDH_end.hour

  wps_nml.write(WPShome+"/"+wps_nml_file, force=True)
  wrf_nml.write(WRFhome+"/"+wrf_nml_file, force=True)
