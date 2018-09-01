import os

def lnsf_wrf(item,WRFwork):
  os.system("ln -sf "+item+" "+WRFwork)

def init_wrf(conf):

  FIXbwrf = conf.get("DEFAULT","FIXbwrf")
  WRFhome = conf.get("DEFAULT","WRFhome")
  WRFwork = conf.get("DEFAULT","WRFwork")  
  WPSwork = conf.get("DEFAULT","WPSwork")

  REALexe = conf.get("exec","real")
  WRFexe = conf.get("exec","wrf")

  namelist = conf.get("wrf","namelist")
  landuse = conf.get("wrf","landuse")
  rrtmg_lw = conf.get("wrf","rrtmg_lw")
  rrtmg_sw = conf.get("wrf","rrtmg_sw")
  vegparm = conf.get("wrf","vegparm")
  soilparm = conf.get("wrf","soilparm")
  genparm = conf.get("wrf","genparm")
  qr_acr_qs = conf.get("wrf","qr_acr_qs")
  qr_acr_qg = conf.get("wrf","qr_acr_qg")
  freezeH2O = conf.get("wrf","freezeH2O")

  WRFwork_exists=os.path.isdir(WRFwork)
  if not WRFwork_exists:
    os.mkdir(WRFwork)

  lnsf_wrf(WRFhome+"/"+namelist,WRFwork)
  lnsf_wrf(WRFhome+"/"+landuse,WRFwork)
  lnsf_wrf(WRFhome+"/"+rrtmg_lw,WRFwork)
  lnsf_wrf(WRFhome+"/"+rrtmg_sw,WRFwork)
  lnsf_wrf(WRFhome+"/"+vegparm,WRFwork)
  lnsf_wrf(WRFhome+"/"+soilparm,WRFwork)
  lnsf_wrf(WRFhome+"/"+genparm,WRFwork)
  lnsf_wrf(WRFhome+"/"+REALexe,WRFwork)
  lnsf_wrf(WRFhome+"/"+WRFexe,WRFwork)
  lnsf_wrf(FIXbwrf+"/"+qr_acr_qs,WRFwork)
  lnsf_wrf(FIXbwrf+"/"+qr_acr_qg,WRFwork)
  lnsf_wrf(FIXbwrf+"/"+freezeH2O,WRFwork)
  lnsf_wrf(WPSwork+"/met_em.*.nc",WRFwork)

def run_real(conf):

  WRFwork = conf.get("DEFAULT","WRFwork")
  REALexe = conf.get("exec","real")

  os.chdir(WRFwork)
  os.system("./"+REALexe+" > real.log")

def run_forecast(conf):

  WRFwork = conf.get("DEFAULT","WRFwork")
  WRFexe = conf.get("exec","wrf")

  os.chdir(WRFwork)
  os.system("./"+WRFexe+" > wrf.log")
