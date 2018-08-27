import os

class geogrid:

  def run_geogrid(self):
    os.system("./geogrid.exe")

class ungrib:

  def link_vtable(self,vtable_path,vtable):
    os.system("ln -s vtable_path+vtable")

  def link_grib(self,grib_path):
    os.system("./link_grib.csh grib_path")

  def run_ungrib(self):
    os.system("./ungrib.exe")

class metgrid:

  def run_metgrid(self):
    os.system("./metgrid.exe")
