import os

class forecast:

  def run_forecast(self):
    os.system("./wrf.exe > wrf.log")
