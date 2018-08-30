import os

def fetch_item(url,item,conf):
  gfs_dir = conf.get("bwrf_data","gfs")
  os.system("wget -O "+gfs_dir+item+" "+url+item)
#  os.system("mv "+item+" "+gfs_dir+item)
