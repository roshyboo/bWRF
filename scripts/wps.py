import os

def lnsf_wps(item,WPSwork):
  os.system("ln -sf "+item+" "+WPSwork)

def init_wps(conf):

  print("-----------------------------------")
  print("------------IN WPS TASK------------")
  print("-----------------------------------")

# Get relevant conf file objects.
  WPShome = conf.get("DEFAULT","WPShome")
  GEOGRIDhome = conf.get("DEFAULT","GEOGRIDhome")
  METGRIDhome = conf.get("DEFAULT","METGRIDhome")
  VTABLEhome = conf.get("wps","VTABLEhome")
  WPSwork = conf.get("DEFAULT","WPSwork")
  GEOGRIDexe = conf.get("exec","geogrid")
  UNGRIBexe = conf.get("exec","ungrib")
  METGRIDexe = conf.get("exec","metgrid")
  gfs = conf.get("bwrf_data","gfs")

  vtable = conf.get("wps","vtable")
  geogrid_table = conf.get("wps","geogrid_table")
  metgrid_table = conf.get("wps","metgrid_table")
  link_script = conf.get("wps","link_script")
  namelist = conf.get("wps","namelist")

# Make the WPS work directory.
  if not os.path.exists(WPSwork):
    os.mkdir(WPSwork)    

# Link files needed for running WPS to WPSwork.
  lnsf_wps(VTABLEhome+"/"+vtable,WPSwork+"/Vtable")
  lnsf_wps(GEOGRIDhome+"/"+geogrid_table,WPSwork+"/GEOGRID.TBL")
  lnsf_wps(METGRIDhome+"/"+metgrid_table,WPSwork+"/METGRID.TBL")
  lnsf_wps(WPShome+"/"+link_script,WPSwork)
  lnsf_wps(WPShome+"/"+GEOGRIDexe,WPSwork)
  lnsf_wps(WPShome+"/"+UNGRIBexe,WPSwork)
  lnsf_wps(WPShome+"/"+METGRIDexe,WPSwork)
  lnsf_wps(WPShome+"/"+namelist,WPSwork)

# Run link script to link the gfs grib files to WPSwork.
  os.chdir(WPSwork)
  os.system("./"+link_script+" "+gfs+"/*grb2*")

def run_geogrid(conf):

# Are geogrid data already available?
# If not, don't run it again.
  FIXbwrf=conf.get("DEFAULT","FIXbwrf")
  WPSwork = conf.get("DEFAULT","WPSwork")
  geogrid_d01_file=conf.get("wps","geogrid_d01_file")

  geogrid_flag = os.path.isfile(FIXbwrf+"/"+geogrid_d01_file)
  if geogrid_flag:
    lnsf_wps(FIXbwrf+"/"+geogrid_d01_file,WPSwork)
  else:
    os.chdir(WPSwork)
    grep_code = os.system("grep 'Successful completion of program geogrid.exe' geogrid.log")
    if grep_code > 0:
      print("Did not see success complete in geogrid.log; will run geogrid.")
      print("-----------------------------------")
      print("----------RUNNING GEOGRID----------")
      print("-----------------------------------")
      GEOGRIDexe = conf.get("exec","geogrid")
      os.system("./"+GEOGRIDexe)

      grep_code = os.system("grep 'Successful completion of program geogrid.exe' geogrid.log")
      if grep_code > 0:
        raise Exception('Program geogrid.exe did not run successfully.')

def run_ungrib(conf):

  WPSwork = conf.get("DEFAULT","WPSwork")
  UNGRIBexe = conf.get("exec","ungrib")

  os.chdir(WPSwork)
  grep_code = os.system("grep 'Successful completion of program ungrib.exe' ungrib.log")
  if grep_code > 0:
    print("Did not see success complete in ungrib.log; will run ungrib.")
    print("-----------------------------------")
    print("----------RUNNING UNGRIB-----------")
    print("-----------------------------------")
    os.system("./"+UNGRIBexe)

    grep_code = os.system("grep 'Successful completion of program ungrib.exe' ungrib.log")
    if grep_code > 0:
      raise Exception('Program ungrib.exe did not run successfully.')

def run_metgrid(conf):

  WPSwork = conf.get("DEFAULT","WPSwork")
  METGRIDexe = conf.get("exec","metgrid")

  os.chdir(WPSwork)
  grep_code = os.system("grep 'Successful completion of program metgrid.exe' metgrid.log")
  if grep_code > 0:
    print("Did not see success complete in metgrid.log; will run metgrid.")
    print("-----------------------------------")
    print("----------RUNNING METGRID----------")
    print("-----------------------------------")
    os.system("./"+METGRIDexe)

    grep_code = os.system("grep 'Successful completion of program metgrid.exe' metgrid.log")
    if grep_code > 0:
      raise Exception('Program metgrid.exe did not run successfully.')
