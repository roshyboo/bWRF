import os

def fetch_item(url,item):
  os.system("wget "+url+item)
