from netCDF4 import Dataset
from datetime import datetime
import os
import numpy as np
import glob
import matplotlib.pyplot as plt
import pandas as pd
from pandas.plotting import register_matplotlib_converters
plt.switch_backend('agg')
from matplotlib.cm import get_cmap
import cartopy.crs as crs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from scipy import interpolate

import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.plots import add_metpy_logo, SkewT
from metpy.units import units

from wrf import CoordPair, vertcross, to_np, getvar, interplevel, smooth2d, get_cartopy, cartopy_xlim, cartopy_ylim, latlon_coords, ALL_TIMES, ll_to_xy

register_matplotlib_converters()

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
    os.mkdir(POSTwork+"/cldfrac")
    os.mkdir(POSTwork+"/low_cldfrac")
    os.mkdir(POSTwork+"/mid_cldfrac")
    os.mkdir(POSTwork+"/high_cldfrac")
    os.mkdir(POSTwork+"/sfc_temp")
    os.mkdir(POSTwork+"/sfcdiags")
    os.mkdir(POSTwork+"/sounding")
    os.mkdir(POSTwork+"/xcdiags")
    os.mkdir(POSTwork+"/xcdiags_rh")
    os.mkdir(POSTwork+"/xcdiags_rh_big")
    os.mkdir(POSTwork+"/700mb")
    os.mkdir(POSTwork+"/500mb")
    os.mkdir(POSTwork+"/300mb")

  lnsf_post(WRFwork+"/wrfout*",POSTwork)

def run_post(conf):

  FIXdir = conf.get("DEFAULT","FIXbwrf")

  sfc_switch = conf.getint("post","plot_sfc")
  cloud_switch = conf.getint("post","plot_clouds")
  sfcdiags_switch = conf.getint("post","plot_sfcdiags")
  xcdiags_switch = conf.getint("post","plot_xcdiags")
  switch_700mb = conf.getint("post","plot_700mb")
  switch_500mb = conf.getint("post","plot_500mb")
  switch_300mb = conf.getint("post","plot_300mb")

#  sfc_switch = 0
#  cloud_switch = 0
#  sfcdiags_switch = 0
#  xcdiags_switch = 0
#  switch_700mb = 0
#  switch_500mb = 0
#  switch_300mb = 0

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

# Reflectivity
  ref = getvar(ncfile, "REFL_10CM", timeidx=ALL_TIMES)
  ref[0,0,0,0]=-21.0 # hack to plot a blank contour plot at the initial time

# Get upper-air quantities.
  p = getvar(ncfile, "pressure", timeidx=ALL_TIMES)
  z = getvar(ncfile, "z", units="dm", timeidx=ALL_TIMES)
  u, v = getvar(ncfile, "uvmet", units="kt", timeidx=ALL_TIMES, meta=False)
  w = getvar(ncfile, "wa", units="m s-1", timeidx=ALL_TIMES)*100. # m/s to cm/s
  rh = getvar(ncfile, "rh", timeidx=ALL_TIMES)
  wspeed = (u**2.0+v**2.0)**0.5
  tc = getvar(ncfile, "tc", timeidx=ALL_TIMES)
  dewT = getvar(ncfile, "td", units="degC", timeidx=ALL_TIMES)
  cloudfrac = getvar(ncfile, "cloudfrac", low_thresh=30., 
                     mid_thresh=955., high_thresh=4500., timeidx=ALL_TIMES)
  total_cloudfrac=np.max(cloudfrac,axis=0)
  low_cloudfrac = cloudfrac[0,:,:,:]
  mid_cloudfrac = cloudfrac[1,:,:,:]
  high_cloudfrac = cloudfrac[2,:,:,:]

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

  div_300 = smooth2d(divergence([u_300*0.514444, v_300*0.514444], dx), 3) # kt to m/s

# Get the sea level pressure
  slp = getvar(ncfile, "slp", timeidx=ALL_TIMES)

# Get the wet bulb temperature
  twb = getvar(ncfile, "twb", units="degC", timeidx=ALL_TIMES)

  slp_levels=np.arange(980,1040,2)
  z_levels=np.arange(504,620,3)
  z_levels_300=np.arange(804,996,6)
  z_levels_700=np.arange(285,351,3)
  t2_levels=np.arange(-20,125,5)
  tc_levels=np.arange(-40,30,2)
  wspeed_levels=np.arange(40,150,10)
  ref_levels=np.arange(-20,60,5)
  rh_levels=np.arange(70,105,5)
  cldfrac_levels=np.arange(0.,1.1,0.1)
  twb_levels=np.arange(0,1,1)

  wup_levels=np.arange(5,55,10)
  wdown_levels=np.arange(-55,-5,10)

  div_levels=np.arange(-110,120,10)

# Get the 10-m u and v wind components.
  u_10m, v_10m = getvar(ncfile, "uvmet10", units="kt", timeidx=ALL_TIMES)
  wind_10m = (u_10m**2.0+v_10m**2.0)**0.5

# Smooth the sea level pressure since it tends to be noisy near the mountains
  smooth_slp = smooth2d(slp, 3)
  twb=twb[:,0,:,:] # lowest model level
  twb=smooth2d(twb,3)

# Get the latitude and longitude points
  lats, lons = latlon_coords(slp)

# Get the cartopy mapping object
  cart_proj = get_cartopy(slp)

  if sfcdiags_switch == 1:

    bx, by = ll_to_xy(ncfile, blat, blon, meta=False, as_int=True)

#   Create a figure
    fig = plt.figure(figsize=(12,9))
    fileout = "precip.png"

  # should be variable ACSNOW
  #  accumulated_snow = getvar(ncfile, "SNOWNC", timeidx=ALL_TIMES)

    liq_equiv = (getvar(ncfile, "RAINC", timeidx=ALL_TIMES) +
                 getvar(ncfile, "RAINNC", timeidx=ALL_TIMES))/25.4 # mm to inches

    liq_equiv_bdu = liq_equiv[:,by,bx]

    prate = np.zeros(numTimes)

    for itime in range(numTimes):

      if itime > 0 and itime < numTimes-1:
        prate[itime] = 0.5*((liq_equiv_bdu[itime+1]+liq_equiv_bdu[itime]) - 
                            (liq_equiv_bdu[itime-1]+liq_equiv_bdu[itime]))
      elif itime == 0:
        prate[0] = liq_equiv_bdu[0] 
      else:
        prate[numTimes-1] = liq_equiv_bdu[numTimes-1]-liq_equiv_bdu[numTimes-2]

    plt.bar(times, prate, width=0.0425, color="black")
    plt.plot(times, liq_equiv_bdu, color="blue")
    plt.xlim(min(times), max(times))
    plt.title("Forecast Precipitation: Boulder, CO")
    plt.xlabel("Time (UTC)")
    plt.ylabel("Running total (line) and rate (bars, inches of liquid)")

    plt.grid(b=True, which="major", axis="both", color="gray", linestyle="--")

    fig.savefig("sfcdiags/"+fileout,bbox_inches='tight')
    plt.close(fig)

#   Create a figure
    fig, ax1 = plt.subplots(figsize=(12,9))
    fileout = "t2m_td2m_wind10m_prate.png"

    t2m = 1.8*(getvar(ncfile, "T2", timeidx=ALL_TIMES)-273.15)+32.
    td2m = getvar(ncfile, "td2", units="degF", timeidx=ALL_TIMES)

    color = 'tab:blue'
    ax1.bar(times, prate, width=0.0425, color="blue")
    ax1.plot(times, liq_equiv_bdu, color="blue")
    ax1.set_ylim(0.,max([max(prate)/0.1+0.01,max(liq_equiv_bdu)+0.1]))
    ax1.set_xlabel("Time (UTC)")
    ax1.set_ylabel("Precipitation liquid amount (in) and rate (bars, in/hr)", color=color)
    ax1.tick_params(axis="y", labelcolor=color)
    ax1.grid(b=True, which="major", axis="x", color="gray", linestyle="--")

    ax2=ax1.twinx()
    ax2.plot(times, wind_10m[:,by,bx], color="black")
    ax2.barbs(times, wind_10m[:,by,bx], u_10m[:,by,bx], v_10m[:,by,bx])

    ax2.plot(times, t2m[:,by,bx], color="red")
    ax2.plot(times, td2m[:,by,bx], color="green")

    ax2.set_xlim(min(times), max(times))
    plt.title("Forecast near-surface variables: Boulder, CO")
    ax2.set_ylabel("2-m T (red) and 2-m Td (green, degF), 10-m wind (black, kt)")

    ax2.grid(b=True, which="major", axis="y", color="gray", linestyle="--")

    fig.tight_layout()
    fig.savefig("sfcdiags/"+fileout,bbox_inches='tight')
    plt.close(fig)

  if xcdiags_switch == 1:

    zinterp = np.arange(550, 880, 10)

    u_xc = vertcross(u, p, levels=zinterp, wrfin=ncfile,
      stagger='u', pivot_point=CoordPair(lat=blat,lon=blon), angle=90., meta=False)

    w_xc = vertcross(w, p, levels=zinterp, wrfin=ncfile,
      stagger='u', pivot_point=CoordPair(lat=blat,lon=blon), angle=90., meta=False)

    tc_xc = vertcross(tc, p, levels=zinterp, wrfin=ncfile,
      stagger='u', pivot_point=CoordPair(lat=blat,lon=blon), angle=90., meta=False)

    rh_xc = vertcross(rh, p, levels=zinterp, wrfin=ncfile,
      stagger='u', pivot_point=CoordPair(lat=blat,lon=blon), angle=90., meta=False)

    nx = np.shape(u_xc)[-1]
    xinterp = np.arange(0,nx,1)

    bx, by = ll_to_xy(ncfile, blat, blon, meta=False, as_int=True)

    for itime in range(numTimes):

#     Create a figure
      fig = plt.figure(figsize=(12,9))
      fileout = "mtnwave_xc"+str(itime).zfill(2)+".png"

      ax=plt.gca()
      ax.set_facecolor("black")

      plt.contourf(xinterp, zinterp, u_xc[itime,:,:], 
        cmap=get_cmap("seismic"), levels=np.arange(-50,55,5))
      plt.colorbar(shrink=.62)

      w_contour = plt.contour(xinterp, zinterp, w_xc[itime,:,:],
        "--", levels=np.arange(-120,-20,20),colors="black")
      plt.clabel(w_contour, inline=1, fontsize=10, fmt="%i")

      w_contour = plt.contour(xinterp, zinterp, w_xc[itime,:,:],
        levels=np.arange(20,120,20),colors="black")
      plt.clabel(w_contour, inline=1, fontsize=10, fmt="%i")

      t_contour = plt.contour(xinterp, zinterp, tc_xc[itime,:,:],
        levels=[-20,-10,0], colors="yellow")
      plt.clabel(t_contour, inline=1, fontsize=10, fmt="%i")

#     Add location of Boulder to plot.
      plt.scatter(bx,np.max(zinterp),c='r',marker='+')

      plt.ylim([np.max(zinterp),np.min(zinterp)])
#      plt.yscale("log")
#      plt.xticks([np.arange(900,475,-25)], ["900", "875", 
#        "850", "825", "800", "775", "750", "725", "700", "675", "650",
#        "625", "600", "575", "550", "525", "500"])

      plt_time=str(dtimes[itime])
      plt_time=plt_time[0:13]

      plt.title(plt_time+"Z fhr "+str(itime).zfill(2)+": zonal wind (fill, kt), temp (yellow lines, degC), VV (black lines, cm/s)")

      fig.savefig("xcdiags/"+fileout,bbox_inches='tight')
      plt.close(fig)

    os.system("convert -delay 90 -dispose background xcdiags/mtnwave_xc*.png -loop 0 xcdiags/mtnwave_xc.gif")

    for itime in range(numTimes):

#     Create a figure
      fig = plt.figure(figsize=(12,9))
      fileout = "mtnwave_xc_rh"+str(itime).zfill(2)+".png"

      ax=plt.gca()
      ax.set_facecolor("black")

      plt.contourf(xinterp, zinterp, rh_xc[itime,:,:],
        cmap=get_cmap("Greens"), levels=rh_levels, extend='both')
      plt.colorbar(shrink=.62)

      skip=2
      plt.quiver(xinterp[::skip], zinterp[::skip], u_xc[itime,::skip,::skip], w_xc[itime,::skip,::skip]/2.,
        scale=500, headwidth=3, color='black', pivot='middle')

      t_contour = plt.contour(xinterp, zinterp, tc_xc[itime,:,:],
        levels=[-20,-10,0], colors="yellow")
      plt.clabel(t_contour, inline=1, fontsize=10, fmt="%i")

#     Add location of Boulder to plot.
      plt.scatter(bx,np.max(zinterp),c='r',marker='+')

      plt.ylim([np.max(zinterp),np.min(zinterp)])

      plt_time=str(dtimes[itime])
      plt_time=plt_time[0:13]

      plt.title(plt_time+"Z fhr "+str(itime).zfill(2)+": rh (fill), temp (yellow lines, degC), wind (arrows)")

      fig.savefig("xcdiags_rh/"+fileout,bbox_inches='tight')
      plt.close(fig)

    os.system("convert -delay 90 -dispose background xcdiags_rh/mtnwave_xc_rh*.png -loop 0 xcdiags_rh/mtnwave_xc_rh.gif")

    zinterp = np.arange(150, 900, 25)

    rh_xc = vertcross(rh, p, levels=zinterp, wrfin=ncfile,
      stagger='u', pivot_point=CoordPair(lat=blat,lon=blon), angle=90., meta=False)

    tc_xc = vertcross(tc, p, levels=zinterp, wrfin=ncfile,
      stagger='u', pivot_point=CoordPair(lat=blat,lon=blon), angle=90., meta=False)

    nx = np.shape(u_xc)[-1]
    xinterp = np.arange(0,nx,1)

    bx, by = ll_to_xy(ncfile, blat, blon, meta=False, as_int=True)

    for itime in range(numTimes):

#     Create a figure
      fig = plt.figure(figsize=(12,9))
      fileout = "mtnwave_xc_rh_big"+str(itime).zfill(2)+".png"

      ax=plt.gca()
      ax.set_facecolor("black")

      plt.contourf(xinterp, zinterp, rh_xc[itime,:,:],
        cmap=get_cmap("Greens"), levels=rh_levels, extend='both')
      plt.colorbar(shrink=.62)

      t_contour = plt.contour(xinterp, zinterp, tc_xc[itime,:,:],
        levels=[-20,-10,0], colors="yellow")
      plt.clabel(t_contour, inline=1, fontsize=10, fmt="%i")

#     Add location of Boulder to plot.
      plt.scatter(bx,np.max(zinterp),c='r',marker='+')

      plt.ylim([np.max(zinterp),np.min(zinterp)])

      plt_time=str(dtimes[itime])
      plt_time=plt_time[0:13]

      plt.title(plt_time+"Z fhr "+str(itime).zfill(2)+": rh (fill), temp (yellow lines, degC)")

      fig.savefig("xcdiags_rh_big/"+fileout,bbox_inches='tight')
      plt.close(fig)

    os.system("convert -delay 90 -dispose background xcdiags_rh_big/mtnwave_xc_rh_big*.png -loop 0 xcdiags_rh_big/mtnwave_xc_rh_big.gif")

  if sfc_switch == 1:

    bx, by = ll_to_xy(ncfile, blat, blon, meta=False, as_int=True)

    for itime in range(numTimes):

      ps = p[itime,:,by,bx]
      ps_temp = ps
      T = tc[itime,:,by,bx]
      Td = dewT[itime,:,by,bx]
      us = u[itime,:,by,bx]
      vs = v[itime,:,by,bx]

      ps = ps[ps_temp >= 100.]
      T = T[ps_temp >= 100.]
      Td = Td[ps_temp >= 100.]
      us = us[ps_temp >= 100.]
      vs = vs[ps_temp >= 100.]

      fig = plt.figure(figsize=(9, 9))
      skew = SkewT(fig, rotation=45)

      fileout = "sounding"+str(itime).zfill(2)+".png"

#     Plot the data using normal plotting functions, in this case using
#     log scaling in Y, as dictated by the typical meteorological plot.
      skew.plot(ps, T, 'r')
      skew.plot(ps, Td, 'g')
      skew.plot_barbs(ps, us, vs)
      skew.ax.set_ylim(1000, 100)
      skew.ax.set_xlim(-40, 60)

#     Calculate LCL height and plot as black dot. Because `p`'s first value is
#     ~1000 mb and its last value is ~250 mb, the `0` index is selected for
#     `p`, `T`, and `Td` to lift the parcel from the surface. If `p` was inverted,
#     i.e. start from low value, 250 mb, to a high value, 1000 mb, the `-1` index
#     should be selected.
      lcl_pressure, lcl_temperature = mpcalc.lcl(ps[0], T[0], Td[0])
      skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black')

#     Calculate full parcel profile and add to plot as black line
      prof = mpcalc.parcel_profile(ps, T[0], Td[0]).to('degC')
      skew.plot(ps, prof, 'k', linewidth=2)

#     Shade areas of CAPE and CIN
#      skew.shade_cin(ps, T, prof, Td)
      [cape, cin] = mpcalc.cape_cin(ps, T, Td, prof)
#      skew.shade_cape(ps, T, prof)

#     An example of a slanted line at constant T -- in this case the 0
#     isotherm
      skew.ax.axvline(0, color='c', linestyle='--', linewidth=2)

#     Add the relevant special lines
      skew.plot_dry_adiabats()
      skew.plot_moist_adiabats()
      skew.plot_mixing_lines()

      plt_time=str(dtimes[itime])
      plt_time=plt_time[0:13]

      plt.title(plt_time+"Z fhr "+str(itime).zfill(2)+": Boulder, CO, CAPE: "+str(round(cape.magnitude))+" J/kg CIN: "+str(round(cin.magnitude))+" J/kg")
      fig.savefig("sounding/"+fileout,bbox_inches='tight')
      plt.close(fig)

    os.system("convert -delay 90 -dispose background sounding/sounding*.png -loop 0 sounding/sounding.gif")

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

    os.system("convert -delay 90 -dispose background sfc/sfc*.png -loop 0 sfc/sfc.gif")

    for itime in range(numTimes):

#   Create a figure
      fig = plt.figure(figsize=(12,9))
      fileout = "sfc_temp_fhr"+str(itime).zfill(2)+".png"

#   Set the GeoAxes to the projection used by WRF
      ax = plt.axes(projection=cart_proj)

#   Add the states and coastlines
      ax.add_feature(states, linewidth=0.8, edgecolor='gray')
      ax.add_feature(COUNTIES, linewidth=0.4, facecolor='none', edgecolor='gray')
      ax.coastlines('50m', linewidth=0.8)

#   2-m air temperature
      t2 = 1.8*(getvar(ncfile, "T2", timeidx=ALL_TIMES) - 273.15)+32.
      plt.contourf(to_np(lons), to_np(lats), to_np(t2[itime,:,:]), transform=crs.PlateCarree(),
             cmap=get_cmap("jet"), levels=t2_levels)

#   Add a color bar
      plt.colorbar(ax=ax, shrink=.62)

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

      plt.title(plt_time+"Z fhr "+str(itime).zfill(2)+": SLP (contours, hPa), 10-m wind (barbs, kt), 2-m T (fill, degF)")
      fig.savefig("sfc_temp/"+fileout,bbox_inches='tight')
      plt.close(fig)

    os.system("convert -delay 90 -dispose background sfc_temp/sfc_temp*.png -loop 0 sfc_temp/sfc_temp.gif")

  if cloud_switch == 1:

    for itime in range(numTimes):

#   Create a figure
      fig = plt.figure(figsize=(12,9))
      fileout = "cldfrac_fhr"+str(itime).zfill(2)+".png"

#   Set the GeoAxes to the projection used by WRF
      ax = plt.axes(projection=cart_proj)

#   Add the states and coastlines
      ax.add_feature(states, linewidth=0.8, edgecolor='gray')
      ax.add_feature(COUNTIES, linewidth=0.4, facecolor='none', edgecolor='gray')
      ax.coastlines('50m', linewidth=0.8)

#   Compute and plot the total cloud fraction.
      plt.contourf(to_np(lons), to_np(lats), to_np(total_cloudfrac[itime,:,:]), transform=crs.PlateCarree(),
             cmap=get_cmap("Greys"), levels=cldfrac_levels, extend='both')

#   Add a color bar
      plt.colorbar(ax=ax, shrink=.62)

#   Add location of Boulder to plot.
      plt.scatter(blon,blat,c='r',marker='+',transform=crs.PlateCarree())

#   Set the map limits.
      ax.set_xlim(cartopy_xlim(cloudfrac))
      ax.set_ylim(cartopy_ylim(cloudfrac))

      plt_time=str(dtimes[itime])
      plt_time=plt_time[0:13]

      plt.title(plt_time+"Z fhr "+str(itime).zfill(2)+": Total cloud fraction (fill)")
      fig.savefig("cldfrac/"+fileout,bbox_inches='tight')
      plt.close(fig)

    os.system("convert -delay 90 -dispose background cldfrac/cldfrac*.png -loop 0 cldfrac/cldfrac.gif")

    for itime in range(numTimes):

#   Create a figure
      fig = plt.figure(figsize=(12,9))
      fileout = "low_cldfrac_fhr"+str(itime).zfill(2)+".png"

#   Set the GeoAxes to the projection used by WRF
      ax = plt.axes(projection=cart_proj)

#   Add the states and coastlines
      ax.add_feature(states, linewidth=0.8, edgecolor='gray')
      ax.add_feature(COUNTIES, linewidth=0.4, facecolor='none', edgecolor='gray')
      ax.coastlines('50m', linewidth=0.8)

#   Compute and plot the cloud fraction.
      plt.contourf(to_np(lons), to_np(lats), to_np(low_cloudfrac[itime,:,:]), transform=crs.PlateCarree(),
             cmap=get_cmap("Greys"), levels=cldfrac_levels, extend='both')

#   Add a color bar
      plt.colorbar(ax=ax, shrink=.62)

#   Add location of Boulder to plot.
      plt.scatter(blon,blat,c='r',marker='+',transform=crs.PlateCarree())

#   Set the map limits.
      ax.set_xlim(cartopy_xlim(cloudfrac))
      ax.set_ylim(cartopy_ylim(cloudfrac))

      plt_time=str(dtimes[itime])
      plt_time=plt_time[0:13]

      plt.title(plt_time+"Z fhr "+str(itime).zfill(2)+": Low cloud fraction (fill)")
      fig.savefig("low_cldfrac/"+fileout,bbox_inches='tight')
      plt.close(fig)

    os.system("convert -delay 90 -dispose background low_cldfrac/low_cldfrac*.png -loop 0 low_cldfrac/low_cldfrac.gif")

    for itime in range(numTimes):

#   Create a figure
      fig = plt.figure(figsize=(12,9))
      fileout = "mid_cldfrac_fhr"+str(itime).zfill(2)+".png"

#   Set the GeoAxes to the projection used by WRF
      ax = plt.axes(projection=cart_proj)

#   Add the states and coastlines
      ax.add_feature(states, linewidth=0.8, edgecolor='gray')
      ax.add_feature(COUNTIES, linewidth=0.4, facecolor='none', edgecolor='gray')
      ax.coastlines('50m', linewidth=0.8)

#   Compute and plot the cloud fraction.
      plt.contourf(to_np(lons), to_np(lats), to_np(mid_cloudfrac[itime,:,:]), transform=crs.PlateCarree(),
             cmap=get_cmap("Greys"), levels=cldfrac_levels, extend='both')

#   Add a color bar
      plt.colorbar(ax=ax, shrink=.62)

#   Add location of Boulder to plot.
      plt.scatter(blon,blat,c='r',marker='+',transform=crs.PlateCarree())

#   Set the map limits.
      ax.set_xlim(cartopy_xlim(cloudfrac))
      ax.set_ylim(cartopy_ylim(cloudfrac))

      plt_time=str(dtimes[itime])
      plt_time=plt_time[0:13]

      plt.title(plt_time+"Z fhr "+str(itime).zfill(2)+": Mid cloud fraction (fill)")
      fig.savefig("mid_cldfrac/"+fileout,bbox_inches='tight')
      plt.close(fig)

    os.system("convert -delay 90 -dispose background mid_cldfrac/mid_cldfrac*.png -loop 0 mid_cldfrac/mid_cldfrac.gif")

    for itime in range(numTimes):

#   Create a figure
      fig = plt.figure(figsize=(12,9))
      fileout = "high_cldfrac_fhr"+str(itime).zfill(2)+".png"

#   Set the GeoAxes to the projection used by WRF
      ax = plt.axes(projection=cart_proj)

#   Add the states and coastlines
      ax.add_feature(states, linewidth=0.8, edgecolor='gray')
      ax.add_feature(COUNTIES, linewidth=0.4, facecolor='none', edgecolor='gray')
      ax.coastlines('50m', linewidth=0.8)

#   Compute and plot the cloud fraction.
      plt.contourf(to_np(lons), to_np(lats), to_np(high_cloudfrac[itime,:,:]), transform=crs.PlateCarree(),
             cmap=get_cmap("Greys"), levels=cldfrac_levels, extend='both')

#   Add a color bar
      plt.colorbar(ax=ax, shrink=.62)

#   Add location of Boulder to plot.
      plt.scatter(blon,blat,c='r',marker='+',transform=crs.PlateCarree())

#   Set the map limits.
      ax.set_xlim(cartopy_xlim(cloudfrac))
      ax.set_ylim(cartopy_ylim(cloudfrac))

      plt_time=str(dtimes[itime])
      plt_time=plt_time[0:13]

      plt.title(plt_time+"Z fhr "+str(itime).zfill(2)+": High cloud fraction (fill)")
      fig.savefig("high_cldfrac/"+fileout,bbox_inches='tight')
      plt.close(fig)

    os.system("convert -delay 90 -dispose background high_cldfrac/high_cldfrac*.png -loop 0 high_cldfrac/high_cldfrac.gif")

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

    os.system("convert -delay 90 -dispose background 500mb/500mb*.png -loop 0 500mb/500mb.gif")

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
             cmap=get_cmap("Greens"), levels=rh_levels, extend='both')

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

    os.system("convert -delay 90 -dispose background 700mb/700mb*.png -loop 0 700mb/700mb.gif")

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

#   300 mb divergence
      plt.contourf(to_np(lons), to_np(lats), to_np(div_300[itime,:,:]), transform=crs.PlateCarree(),
             cmap=get_cmap("seismic"), levels=div_levels)

#   Add a color bar
      plt.colorbar(ax=ax, shrink=.62)

#   Make the 300 mb height contours.
      c_z = plt.contour(to_np(lons), to_np(lats), to_np(z_300[itime,:,:]), transform=crs.PlateCarree(),
             colors="black", levels=z_levels_300)
      plt.clabel(c_z, inline=1, fontsize=10, fmt="%i")

#   Make the 300 mb divergence contours.
#      c_d = plt.contour(to_np(lons), to_np(lats), to_np(div_300[itime,:,:]), transform=crs.PlateCarree(),
#             colors="red", levels=div_levels, linewidths=0.8)
#      plt.clabel(c_d, inline=1, fontsize=10, fontcolor="red", fmt="%i")
#
#      c_c = plt.contour(to_np(lons), to_np(lats), to_np(div_300[itime,:,:]), transform=crs.PlateCarree(),
#             colors="blue", levels=conv_levels, linewidths=0.8)
#      plt.clabel(c_c, inline=1, fontsize=10, fontcolor="blue", fmt="%i")

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

    os.system("convert -delay 90 -dispose background 300mb/300mb*.png -loop 0 300mb/300mb.gif")
