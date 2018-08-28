import os

def fetch_item(url,item,conf):
  data_dir = conf.get("bwrf_data","gfs")
  os.system("wget "+url+item)
  os.system("mv "+item+" "+data_dir+item)
