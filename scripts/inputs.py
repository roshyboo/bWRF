import os

def fetch_item(url,item,input_dir):
  os.system("wget -O "+input_dir+"/"+item+" "+url+"/"+item)

def init_input(conf):
  inputroot=conf.get("bwrf_data","inputroot")
  gfs_input_dir=conf.get("bwrf_data","gfs")

  inputroot_exists=os.path.isdir(inputroot)
  if not inputroot_exists:
    os.mkdir(inputroot)  
    os.mkdir(gfs_input_dir)

def run_input(conf):
# Read the cycle conf.
  gfs_dataset=conf.get("inputs","gfs")
  gfs_grib_item=conf.get("inputs","gfs_gribA")
  gfs_input_dir=conf.get("bwrf_data","gfs")

# Expand the gfs grib item into a list of files.
  run_hours = conf.getfloat("wrf","run_hours")
  interval_hours = conf.getfloat("wrf","interval_seconds")/3600

  ihr = 0
  while ihr <= run_hours:
    fahr = int(ihr)
    conf.set("inputs","fahr",'{:03d}'.format(fahr))
    gfs_grib_item = conf.get("inputs","gfs_gribA")
    fetch_item(gfs_dataset,gfs_grib_item,gfs_input_dir)
    ihr = ihr + interval_hours
