class geogrid

  def run_geogrid(self)
    ./geogrid.exe

class ungrib

  def link_vtable(self,vtable_path,vtable)
    ln -s vtable_path+vtable

  def link_grib(self,grib_path)
    ./link_grib.csh grib_path

  def run_ungrib(self)
    ./ungrib.exe

class metgrid

  def run_metgrid(self)
    ./metgrid.exe
