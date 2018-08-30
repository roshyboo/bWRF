import os

def lnsf_wps(item,WPSwork):
  os.system("ln -sf "+item+" "+WPSwork)

def init_wps(conf):

# Get relevant conf file objects.
  FIXbwrf=conf.get("DEFAULT","FIXbwrf")
  WPShome = conf.get("DEFAULT","WPShome")
  METGRIDhome = conf.get("DEFAULT","METGRIDhome")
  VTABLEhome = conf.get("wps","VTABLEhome")
  WPSwork = conf.get("DEFAULT","WPSwork")
  GEOGRIDexe = conf.get("exec","geogrid")
  UNGRIBexe = conf.get("exec","ungrib")
  METGRIDexe = conf.get("exec","metgrid")
  gfs = conf.get("bwrf_data","gfs")

  geogrid_d01_file=conf.get("wps","geogrid_d01_file")
  vtable = conf.get("wps","vtable")
  metgrid_table = conf.get("wps","metgrid_table")
  link_script = conf.get("wps","link_script")
  namelist = conf.get("wps","namelist")

# Are geogrid data already available?
# If not, don't run it again.
  geogrid_flag = os.path.isfile(FIXbwrf+"/"+geogrid_d01_file)
  if geogrid_flag:
    lnsf_wps(FIXbwrf+"/"+geogrid_d01_file,WPSwork)
  else:
    print("No geogrid data found. You must run geogrid.")

# Link files needed for running WPS to WPSwork.
  lnsf_wps(VTABLEhome+"/"+vtable,WPSwork+"/Vtable")
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

  WPSwork = conf.get("DEFAULT","WPSwork")
  GEOGRIDexe = conf.get("exec","geogrid")

  os.chdir(WPSwork)
  os.system("./"+GEOGRIDexe)

def run_ungrib(conf):

  WPSwork = conf.get("DEFAULT","WPSwork")
  UNGRIBexe = conf.get("exec","ungrib")

  os.chdir(WPSwork)
  os.system("./"+UNGRIBexe)

def run_metgrid(conf):

  WPSwork = conf.get("DEFAULT","WPSwork")
  METGRIDexe = conf.get("exec","metgrid")

  os.chdir(WPSwork)
  os.system("./"+METGRIDexe)
