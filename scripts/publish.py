import subprocess

def run_publish(conf):

  print("-----------------------------------")
  print("---------IN PUBLISH TASK-----------")
  print("-----------------------------------")

  POSTwork = conf.get("DEFAULT","POSTwork")

  files=['/300mb/300mb.gif','/500mb/500mb.gif','/700mb/700mb.gif',
         '/col_rh/col_rh.gif',
         '/sfc/sfc.gif','/sfc_temp/sfc_temp.gif','/sounding/sounding.gif',
         '/sfcdiags/t2m_td2m_wind10m_prate.png','/xcdiags/mtnwave_xc.gif',
         '/xcdiags_rh/mtnwave_xc_rh.gif','/xcdiags_rh_big/mtnwave_xc_rh_big.gif']

  for item in files:
    subprocess.call(["scp", POSTwork+item, "doubloo1@doubloondefender.com:public_html/plots/"])
