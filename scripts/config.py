import ConfigParser

def addVars(conf,var)

  config=ConfigParser.ConfigParser()
  config.read(conf)

  YMDH=config.get("config","YMDH")
  

  config.set("config","aYMDH", YMDH)


