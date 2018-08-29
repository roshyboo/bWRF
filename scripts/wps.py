import os

def run_geogrid(self):
  os.system("./geogrid.exe")

def link_vtable(self,vtable_path,vtable):
  os.system("ln -s vtable_path+vtable")

def link_grib(self,grib_path):
  os.system("./link_grib.csh grib_path")

def run_ungrib(conf):
  WPSbwrf = conf.get("DEFAULT","WPSbwrf")
  UNGRIBbwrf = conf.get("DEFAULT","UNGRIBbwrf")
  vtable = conf.get("wps","vtable")
  os.system("ln -sf "+vtable+" "+)
  os.system("./link_grib.csh grib_path")
  os.system("./ungrib.exe")

def run_metgrid(self):
  os.system("./metgrid.exe")
