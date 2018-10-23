from netCDF4 import Dataset
from datetime import datetime
import os
import numpy as np
import glob
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib.cm import get_cmap
import cartopy.crs as crs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from scipy import interpolate

from wrf import to_np, getvar, interplevel, smooth2d, get_cartopy, cartopy_xlim, cartopy_ylim, latlon_coords, ALL_TIMES

def interpnan(array):

  x = np.arange(0, array.shape[1])
  y = np.arange(0, array.shape[0])

# mask invalid values
  array = np.ma.masked_invalid(array)

  xx, yy = np.meshgrid(x, y)

# get only the valid values
  x1 = xx[~array.mask]
  y1 = yy[~array.mask]

  newarr = array[~array.mask]

  GD1 = interpolate.griddata((x1, y1), newarr.ravel(),
          (xx, yy), method='linear')

  return GD1

def divergence(f,dx):
    """
    Computes the divergence of the vector field f, corresponding to dFx/dx + dFy/dy + ...
    :param f: List of ndarrays, where every item of the list is one dimension of the vector field
    :return: Single ndarray of the same shape as each of the items in f, which corresponds to a scalar field
    """
    num_dims = len(f)
    return np.ufunc.reduce(np.add, [np.gradient(f[i], dx, axis=i)*10.0**5.0 for i in range(num_dims)])

def lnsf_post(item,POSTwork):
  os.system("ln -sf "+item+" "+POSTwork)

def init_post(conf):

  print("-----------------------------------")
  print("-----------IN POST TASK------------")
  print("-----------------------------------")

  WRFwork = conf.get("DEFAULT","WRFwork") 
  POSTwork = conf.get("DEFAULT","POSTwork")

  POSTwork_exists=os.path.isdir(POSTwork)
  if not POSTwork_exists:
    os.mkdir(POSTwork)
    os.mkdir(POSTwork+"/sfc")
    os.mkdir(POSTwork+"/700mb")
    os.mkdir(POSTwork+"/500mb")
    os.mkdir(POSTwork+"/300mb")

  lnsf_post(WRFwork+"/wrfout*",POSTwork)

def run_post(conf):

  FIXdir = conf.get("DEFAULT","FIXbwrf")

  sfc_switch = conf.getint("post","plot_sfc")
  switch_700mb = conf.getint("post","plot_700mb")
  switch_500mb = conf.getint("post","plot_500mb")
  switch_300mb = conf.getint("post","plot_300mb")

  blat = conf.getfloat("post","blat")
  blon = conf.getfloat("post","blon")

  dx = conf.getfloat("wrf","dx")

  POSTwork = conf.get("DEFAULT","POSTwork")
  os.chdir(POSTwork)

  file_wrf = glob.glob("wrfout*")[0]
  print("Processing "+file_wrf)
  init_time = file_wrf[11:21]+" "+file_wrf[22:24]+"Z"

# Download the states and coastlines
  states = cfeature.NaturalEarthFeature(category='cultural', scale='50m', facecolor='none',
                             name='admin_1_states_provinces_shp')
# Get counties.
  reader = shpreader.Reader(FIXdir+'/shapefiles/countyp010g.shp')
  counties = list(reader.geometries())
  COUNTIES = cfeature.ShapelyFeature(counties, crs.PlateCarree())

# Open the NetCDF file
  ncfile = Dataset(file_wrf)

# Get the times and convert to datetimes
  times = getvar(ncfile, "Times", timeidx=ALL_TIMES, meta=False)
  dtimes = [datetime.strptime(str(atime), '%Y-%m-%dT%H:%M:%S.000000000') for atime in times]
  numTimes = len(times)

# Map factors
# mfmx = getvar(ncfile, "MAPFAC_MX", timeidx=ALL_TIMES)
# mfmy = getvar(ncfile, "MAPFAC_MY", timeidx=ALL_TIMES)

# Reflectivity
  ref = getvar(ncfile, "REFL_10CM", timeidx=ALL_TIMES)
  ref[0,0,0,0]=-19.0 # hack to plot a blank contour plot at the initial time

# Get upper-air quantities.
  p = getvar(ncfile, "pressure", timeidx=ALL_TIMES)
  z = getvar(ncfile, "z", units="dm", timeidx=ALL_TIMES)
  u, v = getvar(ncfile, "uvmet", units="kt", timeidx=ALL_TIMES, meta=False)
  w = getvar(ncfile, "wa", units="m s-1", timeidx=ALL_TIMES)*100. # m/s to cm/s
  rh = getvar(ncfile, "rh", timeidx=ALL_TIMES)
  wspeed = (u**2.0+v**2.0)**0.5
  tc = getvar(ncfile, "tc", timeidx=ALL_TIMES)

# Interpolate
  z_500 = smooth2d(interplevel(z, p, 500), 3)
  tc_500 = smooth2d(interplevel(tc, p, 500), 3)
  u_500 = interplevel(u, p, 500)
  v_500 = interplevel(v, p, 500)
  wspeed_500 = interplevel(wspeed, p, 500)

  z_300 = smooth2d(interplevel(z, p, 300), 3)
  u_300 = interplevel(u, p, 300)
  v_300 = interplevel(v, p, 300)

  z_700 = interplevel(z, p, 700)
  tc_700 = interplevel(tc, p, 700)
  u_700 = interplevel(u, p, 700)
  v_700 = interplevel(v, p, 700)
  w_700 = interplevel(w, p, 700)
  rh_700 = interplevel(rh, p, 700)

  if switch_700mb == 1:

# Interpolate over NaNs.

    for itime in range(numTimes):

      z_700[itime,:,:] = interpnan(z_700[itime,:,:])
      tc_700[itime,:,:] = interpnan(tc_700[itime,:,:])
      rh_700[itime,:,:] = interpnan(rh_700[itime,:,:])
      u_700[itime,:,:] = interpnan(u_700[itime,:,:])
      v_700[itime,:,:] = interpnan(v_700[itime,:,:])
      w_700[itime,:,:] = interpnan(w_700[itime,:,:])

    z_700 = smooth2d(z_700, 3)
    tc_700 = smooth2d(tc_700, 3)
    w_700 = smooth2d(w_700, 3)

# dx=30000.*mfmx # mfmx = map factor on mass grid in x direction
# dy=30000.*mfmy # mfmy = map factor on mass grid in y direction
  div_300 = smooth2d(divergence([u_300*0.514444, v_300*0.514444], dx), 3) # kt to m/s

# Get the sea level pressure
  slp = getvar(ncfile, "slp", timeidx=ALL_TIMES)

# Get the wet bulb temperature
  twb = getvar(ncfile, "twb", units="degC", timeidx=ALL_TIMES)

  slp_levels=np.arange(980,1040,2)
  z_levels=np.arange(504,620,3)
  z_levels_300=np.arange(804,996,6)
  z_levels_700=np.arange(285,351,3)
  tc_levels=np.arange(-40,30,2)
  wspeed_levels=np.arange(40,150,10)
  ref_levels=np.arange(-20,60,5)
  rh_levels=np.arange(70,105,5)
  twb_levels=np.arange(0,1,1)

  wup_levels=np.arange(1,50,2)
  wdown_levels=np.arange(-49,0,2)

  div_levels=np.arange(2,20,2)
  conv_levels=np.arange(-20,0,2)

# Get the 10-m u and v wind components.
  u_10m, v_10m = getvar(ncfile, "uvmet10", units="kt", timeidx=ALL_TIMES)

# Smooth the sea level pressure since it tends to be noisy near the mountains
  smooth_slp = smooth2d(slp, 3)
  twb=twb[:,0,:,:] # lowest model level
  twb=smooth2d(twb,3)

# Get the latitude and longitude points
  lats, lons = latlon_coords(slp)

# Get the cartopy mapping object
  cart_proj = get_cartopy(slp)

  if sfc_switch == 1:

    for itime in range(numTimes):

#   Create a figure
      fig = plt.figure(figsize=(12,9))
      fileout = "sfc_fhr"+str(itime).zfill(2)+".png"

#   Set the GeoAxes to the projection used by WRF
      ax = plt.axes(projection=cart_proj)

#   Add the states and coastlines
      ax.add_feature(states, linewidth=0.8, edgecolor='gray')
      ax.add_feature(COUNTIES, linewidth=0.4, facecolor='none', edgecolor='gray')
      ax.coastlines('50m', linewidth=0.8)

#   Reflectivity at the lowest model level.
      plt.contourf(to_np(lons), to_np(lats), to_np(ref[itime,0,:,:]), transform=crs.PlateCarree(),
             cmap=get_cmap("jet"), levels=ref_levels)

#   Add a color bar
      plt.colorbar(ax=ax, shrink=.62)

#   Contour the wetbulb temperature at 0 degC.
      c_p = plt.contour(to_np(lons), to_np(lats), to_np(twb[itime,:,:]), transform=crs.PlateCarree(),
             colors="red", levels=twb_levels)
      plt.clabel(c_p, inline=1, fontsize=10, fmt="%i")

#   Make the contour outlines and filled contours for the smoothed sea level pressure.
      c_p = plt.contour(to_np(lons), to_np(lats), to_np(smooth_slp[itime,:,:]), transform=crs.PlateCarree(),
             colors="black", levels=slp_levels)
      plt.clabel(c_p, inline=1, fontsize=10, fmt="%i")

#   Add location of Boulder to plot.
      plt.scatter(blon,blat,c='r',marker='+',transform=crs.PlateCarree())

#   Add the 10-m wind barbs, only plotting every other data point.
      skip=2
      plt.barbs(to_np(lons[::skip,::skip]), to_np(lats[::skip,::skip]), to_np(u_10m[itime, ::skip, ::skip]),
          to_np(v_10m[itime, ::skip, ::skip]), transform=crs.PlateCarree(), length=5.25, linewidth=0.5)

#   Set the map limits.
      ax.set_xlim(cartopy_xlim(smooth_slp))
      ax.set_ylim(cartopy_ylim(smooth_slp))

      plt_time=str(dtimes[itime])
      plt_time=plt_time[0:13]

      plt.title(plt_time+"Z fhr "+str(itime).zfill(2)+": SLP (fill, hPa), 10-m wind (barbs, kt), LML Ref (fill, dBZ), LML WBT (red, degC)")
      fig.savefig("sfc/"+fileout,bbox_inches='tight')
      plt.close(fig)

    os.system("convert -delay 25 -dispose background sfc/sfc*.png -loop 0 sfc/sfc.gif")

  if switch_500mb == 1:

    for itime in range(numTimes):

#   Create a figure
      fig = plt.figure(figsize=(12,9))
      fileout = "500mb_fhr"+str(itime).zfill(2)+".png"

#   Set the GeoAxes to the projection used by WRF
      ax = plt.axes(projection=cart_proj)

#   Add the states and coastlines
      ax.add_feature(states, linewidth=0.8, edgecolor='gray')
      ax.add_feature(COUNTIES, linewidth=0.4, facecolor='none', edgecolor='gray')
      ax.coastlines('50m', linewidth=0.8)

#   wind color fill
      if np.max(wspeed_500[itime,:,:]) > np.min(wspeed_levels):
        plt.contourf(to_np(lons), to_np(lats), to_np(wspeed_500[itime,:,:]), transform=crs.PlateCarree(),
             cmap=get_cmap("rainbow"), levels=wspeed_levels)

#   Make the 500 mb height contours.
      c_z = plt.contour(to_np(lons), to_np(lats), to_np(z_500[itime,:,:]), transform=crs.PlateCarree(),
             colors="black", levels=z_levels)
      plt.clabel(c_z, inline=1, fontsize=10, fmt="%i")

#   Make the 500 mb temp contours.
      c_t = plt.contour(to_np(lons), to_np(lats), to_np(tc_500[itime,:,:]), transform=crs.PlateCarree(),
             colors="red", levels=tc_levels)
      plt.clabel(c_t, inline=1, fontsize=10, fontcolor="red", fmt="%i")

#   Add a color bar
#   plt.colorbar(ax=ax, shrink=.62)

#   Add location of Boulder to plot.
      plt.scatter(blon,blat,c='r',marker='+',transform=crs.PlateCarree())

#   Add the 500mb wind barbs, only plotting every other data point.
      skip=3
      plt.barbs(to_np(lons[::skip,::skip]), to_np(lats[::skip,::skip]), to_np(u_500[itime, ::skip, ::skip]),
          to_np(v_500[itime, ::skip, ::skip]), transform=crs.PlateCarree(), length=5.25, linewidth=0.5)

#   Set the map limits.
      ax.set_xlim(cartopy_xlim(z_500))
      ax.set_ylim(cartopy_ylim(z_500))

      plt_time=str(dtimes[itime])
      plt_time=plt_time[0:13]

      plt.title(plt_time+"Z fhr "+str(itime).zfill(2)+": 500-mb height (black, dm), temp (red, degC), and wind (fill/barbs, kt)")
      fig.savefig("500mb/"+fileout,bbox_inches='tight')
      plt.close(fig)

    os.system("convert -delay 25 -dispose background 500mb/500mb*.png -loop 0 500mb/500mb.gif")

  if switch_700mb == 1:

    for itime in range(numTimes):

#   Create a figure
      fig = plt.figure(figsize=(12,9))
      fileout = "700mb_fhr"+str(itime).zfill(2)+".png"

#   Set the GeoAxes to the projection used by WRF
      ax = plt.axes(projection=cart_proj)

#   Add the states and coastlines
      ax.add_feature(states, linewidth=0.8, edgecolor='gray')
      ax.add_feature(COUNTIES, linewidth=0.4, facecolor='none', edgecolor='gray')
      ax.coastlines('50m', linewidth=0.8)

#   rh color fill
      plt.contourf(to_np(lons), to_np(lats), to_np(rh_700[itime,:,:]), transform=crs.PlateCarree(),
             cmap=get_cmap("Greens"), levels=rh_levels)

#   Add a color bar
      plt.colorbar(ax=ax, shrink=.62)

#   Make the 700 mb height contours.
      c_z = plt.contour(to_np(lons), to_np(lats), to_np(z_700[itime,:,:]), transform=crs.PlateCarree(),
             colors="black", levels=z_levels_700)
      plt.clabel(c_z, inline=1, fontsize=10, fmt="%i")

#   Make the 700 mb temp contours.
      c_t = plt.contour(to_np(lons), to_np(lats), to_np(tc_700[itime,:,:]), transform=crs.PlateCarree(),
             colors="red", levels=tc_levels)
      plt.clabel(c_t, inline=1, fontsize=10, fontcolor="red", fmt="%i")

#   Make the 700 mb VV contours.
      c_d = plt.contour(to_np(lons), to_np(lats), to_np(w_700[itime,:,:]), transform=crs.PlateCarree(),
             colors="magenta", levels=wup_levels, linewidths=0.9)
      plt.clabel(c_d, inline=1, fontsize=10, fontcolor="magenta", fmt="%i")

      c_c = plt.contour(to_np(lons), to_np(lats), to_np(w_700[itime,:,:]), transform=crs.PlateCarree(),
             colors="blue", levels=wdown_levels, linewidths=0.9)
      plt.clabel(c_c, inline=1, fontsize=10, fontcolor="blue", fmt="%i")

#   Add location of Boulder to plot.
      plt.scatter(blon,blat,c='r',marker='+',transform=crs.PlateCarree())

#   Add the 700mb wind barbs, only plotting every other data point.
      skip=3
      plt.barbs(to_np(lons[::skip,::skip]), to_np(lats[::skip,::skip]), to_np(u_700[itime, ::skip, ::skip]),
          to_np(v_700[itime, ::skip, ::skip]), transform=crs.PlateCarree(), length=5.25, linewidth=0.5)

#   Set the map limits.
      ax.set_xlim(cartopy_xlim(z_700))
      ax.set_ylim(cartopy_ylim(z_700))

      plt_time=str(dtimes[itime])
      plt_time=plt_time[0:13]

      plt.title(plt_time+"Z fhr "+str(itime).zfill(2)+": 700-mb hgt (black, dm), T (red, degC), wind (barbs, kt), VV (cm/s), rh (fill, %)")
      fig.savefig("700mb/"+fileout,bbox_inches='tight')
      plt.close(fig)

    os.system("convert -delay 25 -dispose background 700mb/700mb*.png -loop 0 700mb/700mb.gif")

  if switch_300mb == 1:

    for itime in range(numTimes):

#   Create a figure
      fig = plt.figure(figsize=(12,9))
      fileout = "300mb_fhr"+str(itime).zfill(2)+".png"

#   Set the GeoAxes to the projection used by WRF
      ax = plt.axes(projection=cart_proj)

#   Add the states and coastlines
      ax.add_feature(states, linewidth=0.8, edgecolor='gray')
      ax.add_feature(COUNTIES, linewidth=0.4, facecolor='none', edgecolor='gray')
      ax.coastlines('50m', linewidth=0.8)

#   Make the 300 mb height contours.
      c_z = plt.contour(to_np(lons), to_np(lats), to_np(z_300[itime,:,:]), transform=crs.PlateCarree(),
             colors="black", levels=z_levels_300)
      plt.clabel(c_z, inline=1, fontsize=10, fmt="%i")

#   Make the 300 mb divergence contours.
      c_d = plt.contour(to_np(lons), to_np(lats), to_np(div_300[itime,:,:]), transform=crs.PlateCarree(),
             colors="red", levels=div_levels, linewidths=0.8)
      plt.clabel(c_d, inline=1, fontsize=10, fontcolor="red", fmt="%i")

      c_c = plt.contour(to_np(lons), to_np(lats), to_np(div_300[itime,:,:]), transform=crs.PlateCarree(),
             colors="blue", levels=conv_levels, linewidths=0.8)
      plt.clabel(c_c, inline=1, fontsize=10, fontcolor="blue", fmt="%i")

#   Add location of Boulder to plot.
      plt.scatter(blon,blat,c='r',marker='+',transform=crs.PlateCarree())

#   Add the 300mb wind barbs, only plotting every other data point.
      skip=3
      plt.barbs(to_np(lons[::skip,::skip]), to_np(lats[::skip,::skip]), to_np(u_300[itime, ::skip, ::skip]),
          to_np(v_300[itime, ::skip, ::skip]), transform=crs.PlateCarree(), length=5.25, linewidth=0.5)

#   Set the map limits.
      ax.set_xlim(cartopy_xlim(z_300))
      ax.set_ylim(cartopy_ylim(z_300))

      plt_time=str(dtimes[itime])
      plt_time=plt_time[0:13]

      plt.title(plt_time+"Z fhr "+str(itime).zfill(2)+": 300-mb height (black, dm), wind (barbs, kt), and divergence x 10^5 (red/blue, s^-1)")
      fig.savefig("300mb/"+fileout,bbox_inches='tight')
      plt.close(fig)

    os.system("convert -delay 25 -dispose background 300mb/300mb*.png -loop 0 300mb/300mb.gif")
