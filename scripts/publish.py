import subprocess, time

def run_publish(conf):

  print("-----------------------------------")
  print("---------IN PUBLISH TASK-----------")
  print("-----------------------------------")

  POSTwork = conf.get("DEFAULT","POSTwork")

  files=['/300mb/300mb.gif','/500mb/500mb.gif','/700mb/700mb.gif',
         '/cldfrac/cldfrac.gif','/low_cldfrac/low_cldfrac.gif',
         '/mid_cldfrac/mid_cldfrac.gif','/high_cldfrac/high_cldfrac.gif',
         '/sfc/sfc.gif','/sfc_temp/sfc_temp.gif','/sounding/sounding.gif',
         '/sfcdiags/t2m_td2m_wind10m_prate.png','/xcdiags/mtnwave_xc.gif',
         '/xcdiags_rh/mtnwave_xc_rh.gif','/xcdiags_rh_big/mtnwave_xc_rh_big.gif']

  for ix,item in enumerate(files):
    if ix == 10:
      time.sleep(30)
    subprocess.call(["scp", POSTwork+item, "doubloo1@doubloondefender.com:public_html/plots/"])
