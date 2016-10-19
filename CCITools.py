from __future__ import division
from analyseCCI import CCI
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial
import sys
import numpy.ma as ma
import math
from netCDF4 import Dataset
import colorsys
import time

def writeCCI(path, data, targetGrid, primary, NOAA = True):
        
    ncOut = Dataset(path, "w", format="NETCDF4")

    along_track  = ncOut.createDimension("along_track" , targetGrid.lon.shape[0])
    across_track = ncOut.createDimension("across_track", targetGrid.lon.shape[1])

    lat = ncOut.createVariable("lat","f8",("along_track", "across_track"))
    lat[:,:] = targetGrid.lat

    lon = ncOut.createVariable("lon","f8",("along_track", "across_track"))
    lon[:,:] = targetGrid.lon

    if primary:
        times = ncOut.createVariable("time","f8",("along_track", "across_track"))
        times[:,:] = data[:,:,0]

        solZen = ncOut.createVariable("solar_zenith_view_no1","f8",("along_track", "across_track"))
        solZen[:,:] = data[:,:,1]

        satZen = ncOut.createVariable("satellite_zenith_view_no1","f8",("along_track", "across_track"))
        satZen[:,:] = data[:,:,2]

        cot = ncOut.createVariable("cot","f8",("along_track", "across_track"))
        cot[:,:] = data[:,:,3]

        cer = ncOut.createVariable("cer","f8",("along_track", "across_track"))
        cer[:,:] = data[:,:,4]
        
        cwp = ncOut.createVariable("cwp","f8",("along_track", "across_track"))
        cwp[:,:] = data[:,:,5]
        
        ctp = ncOut.createVariable("ctp","f8",("along_track", "across_track"))
        ctp[:,:] = data[:,:,6]
        
        ctp_corrected = ncOut.createVariable("ctp_corrected","f8",("along_track", "across_track"))
        ctp_corrected[:,:] = data[:,:,7]
        
        cth = ncOut.createVariable("cth","f8",("along_track", "across_track"))
        cth[:,:] = data[:,:,8]
        
        cth_corrected = ncOut.createVariable("cth_corrected","f8",("along_track", "across_track"))
        cth_corrected[:,:] = data[:,:,9]
        
        ctt = ncOut.createVariable("ctt","f8",("along_track", "across_track"))
        ctt[:,:] = data[:,:,10]

        ctt_corrected = ncOut.createVariable("ctt_corrected","f8",("along_track", "across_track"))
        ctt_corrected[:,:] = data[:,:,11]

        cc_total = ncOut.createVariable("cc_total","f8",("along_track", "across_track"))
        cc_total[:,:] = data[:,:,12]
        
        cccot = ncOut.createVariable("cccot_pre","f8",("along_track", "across_track"))
        cccot[:,:] = data[:,:,13]
        
        phase = ncOut.createVariable("phase","f8",("along_track", "across_track"))
        phase[:,:] = data[:,:,14]

    else:

        if NOAA:
            albedo1 = ncOut.createVariable("albedo_in_channel_no_1","f8",("along_track", "across_track"))
            albedo1[:,:] = data[:,:,0]
            
            albedo2 = ncOut.createVariable("albedo_in_channel_no_2","f8",("along_track", "across_track"))
            albedo2[:,:] = data[:,:,1]

            albedo3 = ncOut.createVariable("albedo_in_channel_no_3","f8",("along_track", "across_track"))
            albedo3[:,:] = data[:,:,2]

            reflectance1 = ncOut.createVariable("reflectance_in_channel_no_1","f8",("along_track", "across_track"))
            reflectance1[:,:] = data[:,:,3]

            reflectance2 = ncOut.createVariable("reflectance_in_channel_no_2","f8",("along_track", "across_track"))
            reflectance2[:,:] = data[:,:,4]

            reflectance3 = ncOut.createVariable("reflectance_in_channel_no_3","f8",("along_track", "across_track"))
            reflectance3[:,:] = data[:,:,5]

            bt4 = ncOut.createVariable("brightness_temperature_in_channel_no_4","f8",("along_track", "across_track"))
            bt4[:,:] = data[:,:,6]

            bt5 = ncOut.createVariable("brightness_temperature_in_channel_no_5","f8",("along_track", "across_track"))
            bt5[:,:] = data[:,:,7]
        
            bt6 = ncOut.createVariable("brightness_temperature_in_channel_no_6","f8",("along_track", "across_track"))
            bt6[:,:] = data[:,:,8]

        else:

            albedo1 = ncOut.createVariable("albedo_in_channel_no_1","f8",("along_track", "across_track"))
            albedo1[:,:] = data[:,:,0]
            
            albedo2 = ncOut.createVariable("albedo_in_channel_no_2","f8",("along_track", "across_track"))
            albedo2[:,:] = data[:,:,1]

            albedo6 = ncOut.createVariable("albedo_in_channel_no_6","f8",("along_track", "across_track"))
            albedo6[:,:] = data[:,:,2]

            reflectance1 = ncOut.createVariable("reflectance_in_channel_no_1","f8",("along_track", "across_track"))
            reflectance1[:,:] = data[:,:,3]

            reflectance2 = ncOut.createVariable("reflectance_in_channel_no_2","f8",("along_track", "across_track"))
            reflectance2[:,:] = data[:,:,4]

            reflectance6 = ncOut.createVariable("reflectance_in_channel_no_6","f8",("along_track", "across_track"))
            reflectance6[:,:] = data[:,:,5]

            bt20 = ncOut.createVariable("brightness_temperature_in_channel_no_20","f8",("along_track", "across_track"))
            bt20[:,:] = data[:,:,6]

            bt31 = ncOut.createVariable("brightness_temperature_in_channel_no_31","f8",("along_track", "across_track"))
            bt31[:,:] = data[:,:,7]
        
            bt32 = ncOut.createVariable("brightness_temperature_in_channel_no_32","f8",("along_track", "across_track"))
            bt32[:,:] = data[:,:,8]

    ncOut.close()

def resampleCCI(sourceCCI, targetGrid, lat_in = None, lon_in = None, NOAA = True):
    # aim: resample CCI L2 orbit data from orbit projection to regular lat/lon grid
    # grid dimension and resolution is user defined
    # approach: average all pixel values whose centre are within grid box

    if lat_in is not None:
        primary = False
    else:
        primary = True

    targetPoints = np.asarray(zip(targetGrid.lon.ravel(), targetGrid.lat.ravel()))
    if primary:
        sourcePoints = np.asarray(zip(sourceCCI.lon.ravel(), sourceCCI.lat.ravel()))
    else:
        sourcePoints = np.asarray(zip(lon_in.ravel(), lat_in.ravel()))
    tree = spatial.cKDTree(targetPoints)

    # variables to be resampled
    if primary:
        variables = ["time", "solar_zenith_view_no1", "satellite_zenith_view_no1", "cot", "cer", "cwp",
                     "ctp", "ctp_corrected", "cth", "cth_corrected", "ctt", "ctt_corrected", "cc_total",
                     "cccot_pre", "phase"]
    else:
        if NOAA:
            variables = ["albedo_in_channel_no_1", "albedo_in_channel_no_2", "albedo_in_channel_no_3",
                         "reflectance_in_channel_no_1", "reflectance_in_channel_no_2", "reflectance_in_channel_no_3",
                         "brightness_temperature_in_channel_no_4", "brightness_temperature_in_channel_no_5",
                         "brightness_temperature_in_channel_no_6"]    
        else:
            variables = ["albedo_in_channel_no_1", "albedo_in_channel_no_2", "albedo_in_channel_no_6",
                         "reflectance_in_channel_no_1", "reflectance_in_channel_no_2", "reflectance_in_channel_no_6",
                         "brightness_temperature_in_channel_no_20", "brightness_temperature_in_channel_no_31",
                         "brightness_temperature_in_channel_no_32"]    

    if primary:
        sourceValues = np.zeros((np.prod(sourceCCI.time.shape), len(variables)))
    else:
        sourceValues = np.zeros((np.prod(sourceCCI.cot_ap.shape), len(variables)))

    for var_i, var in enumerate(variables):
        sourceValues[:, var_i] = getattr(sourceCCI, var).ravel()
                 
    returnValues = np.zeros((targetPoints.shape[0], len(variables)))
    returnCount  = np.zeros((targetPoints.shape[0], len(variables)))

    i = 0
    for lon, lat in sourcePoints:
        nn = tree.query((lon, lat), k = 1)
        targetIndex = nn[1]
        for var_i, var in enumerate(variables):
            if sourceValues[i, var_i] > -32767.0:
                returnCount [targetIndex, var_i] += 1
                returnValues[targetIndex, var_i] += sourceValues[i, var_i] 
        i += 1
        if i % 250000 == 0:
            print ('{0}\r'.format(round(i / sourcePoints.shape[0] * 100, 1)))   

    returnCountMasked = ma.masked_values(returnCount, 0.)
    out = returnValues / returnCountMasked
    out = out.reshape(targetGrid.lat.shape[0], targetGrid.lat.shape[1], len(variables))

    return out

def minMax(a):
    try:
        return min(a), max(a)
    except:
        return a.min(), a.max()

def plotVariable(CCIpri, CCIsec,
                 lat_0, lon_0,
                 variable = "cot", input = None, plotInput = False,
                 width = 5000000, height = 5000000,
                 res = 'l',
                 llcrnrlat = 50, urcrnrlat = 80, llcrnrlon = -180, urcrnrlon = -150,
                 cmin = None, cmax = None,
                 mask = None):    
    
    # define map projection
    map = Basemap(width = width, height = height,
                  resolution='l', projection = 'stere',
                  lat_ts = lat_0, lat_0 = lat_0, lon_0 = lon_0)

    # define figure size
    fig1 = plt.figure(figsize = (10, 10))
    # define # subplots
    ax = fig1.add_subplot(111)

    # draw coasts and fill continents
    map.drawcoastlines(linewidth = 0.5)
    fillContinents = map.fillcontinents(color='#C0C0C0', lake_color='#7093DB', zorder = 0)

    # choose secondary dataset for reflectances/BTs, else take primary
    if plotInput:
        var = input
    else:
        if "reflectance" in variable or "brightness" in variable:
            var = getattr(CCIsec, variable)
        else:
            var = getattr(CCIpri, variable)

    if mask is not None:
        print "    masking"
        var.mask = mask 

    # plot the data, adding a colorbar
    forColorbar = map.pcolormesh(CCIpri.lon, CCIpri.lat, var, latlon = True)
    if cmin is not None and cmax is not None:
        forColorbar.set_clim(vmin = cmin, vmax = cmax)
    map.colorbar(forColorbar, location = 'right')
    

    # add centre point
    x, y = map(lon_0, lat_0)
    map.plot(x, y, 'kx', markersize=15)

    # draw grid lines
    gridSpacing = 5.
    map.drawparallels(
        np.arange(llcrnrlat, urcrnrlat, gridSpacing),
        color = 'black', linewidth = 0.5,
        labels=[True, False, False, False])
    map.drawmeridians(
        np.arange(-180, urcrnrlon, gridSpacing),
        color = '0.25', linewidth = 0.5,
        labels=[False, False, False, True])

    # alternative map projection
    # map = Basemap(projection = 'stere', lat_0 = lat_0, lon_0 = lon_0, 
    #               llcrnrlat = llcrnrlat, urcrnrlat = urcrnrlat, llcrnrlon = llcrnrlon, urcrnrlon = urcrnrlon,
    #               resolution = 'l'
    #               )

def plotCCI(priN18, priMYD, boundingBox, centrePoint, variable, platform,
            secN18 = None, secMYD = None, input = None, colourMin = None, 
            colourMax = None, mask = None):

    # projection type: stereographic
    lat_0     = centrePoint[0]
    lon_0     = centrePoint[1]
    llcrnrlat = boundingBox[2] 
    urcrnrlat = boundingBox[3] 
    llcrnrlon = boundingBox[0] 
    urcrnrlon = boundingBox[1]         

    if input is None:
        plotInput = False
    else:
        plotInput = True

    if colourMin is not None and colourMax is not None:
        cmin = colourMin
        cmax = colourMax
    elif colourMin is not None and colourMax is None:
        cmin = colourMin
        if plotInput:
            cmax = input.max()
        else:
            cmax = variable.max()
    elif colourMin is None and colourMax is not None:
        if plotInput:
            cmin = input.min()
        else:
            cmin = variable.min()
        cmax = colourMax
    else:
        cmin = None
        cmax = None
        
    if platform is 'NOAA18':
        primary = priN18
        secondary = secN18
    elif platform is 'MYD':
        primary = priMYD
        secondary = secMYD
    elif platform is 'MYDResampled':
        primary = MYDResampled
        secondary = secMYD # resampled secondary data not yet available
    elif platform is 'MYDResampledMinusNOAA18':
        primary = priN18
        secondary = secN18
    else:
        print("Define satellite platform to plot: NOAA18, MYD, or MYDResampled. Given: " + platform + ".")
        sys.exit()

    print("Plotting " + variable + " for " + platform + ".")
    plotVariable(primary, secondary,
                 lat_0, lon_0,
                 variable = variable, input = input, plotInput = plotInput,
                 width = 4000000, height = 4000000,
                 res = 'l',
                 llcrnrlat = llcrnrlat, urcrnrlat = urcrnrlat, llcrnrlon = -180, urcrnrlon = urcrnrlon,
                 cmin = cmin, cmax = cmin,
                 mask = mask)
    plt.draw()

def buildRGB(primaryData, secondaryData, platform):

    # set channel numbers according to input platform
    if platform is "MYD":
        ch4 = "brightness_temperature_in_channel_no_20"
        ch5 = "brightness_temperature_in_channel_no_31"
    elif platform is "N18":
        ch4 = "brightness_temperature_in_channel_no_4"
        ch5 = "brightness_temperature_in_channel_no_5"
    else:
        print "Error: select one of [MYD, N18] as input platform."
        sys.exit()
        
    # get variables if attributes not available
    if not hasattr(primaryData, "solar_zenith_view_no1"):
        primaryData.getAllVariables()
    if not hasattr(secondaryData, "albedo_in_channel_no_1"):
        secondaryData.getAllVariables()        

    # define new numpy cosine function that accepts degrees as argument
    np.cosd = lambda x : np.cos( np.deg2rad(x) )
        
    # get data and mask negative values
    sunzen = primaryData.solar_zenith_view_no1
    ir = getattr(secondaryData, ch5)  
    ir = ir.view(ma.MaskedArray) # explicitly add mask to variable, as no missing data were found
    ir.mask = ma.masked_less(getattr(secondaryData, ch5), 0.).mask
    red  = getattr(secondaryData, ch4)
    red.mask = ma.masked_less(getattr(secondaryData, ch4), 0.).mask
    # extract solar component from mixed 3.7 um channel
    red = ((red - ir) / np.cosd(sunzen) + 20. + ((ir - 302.).clip(min = 0.) / 30. * 40.)).clip(min = 0.)
    green = secondaryData.reflectance_in_channel_no_2
    green.mask = ma.masked_less(secondaryData.reflectance_in_channel_no_2, 0.).mask
    green /= np.cosd(sunzen)
    blue   = secondaryData.reflectance_in_channel_no_1
    blue.mask = ma.masked_less(secondaryData.reflectance_in_channel_no_1, 0.).mask
    blue /= np.cosd(sunzen)

    # scale into bit space - legacy from IDL code
    red   = ((red   / 150.      ).clip(min = 0.)).clip(max = 1.) * 255.
    green = ((green / 138. * 0.8).clip(min = 0.)).clip(max = 1.) * 255.
    blue  = ((blue  / 122. * 0.5).clip(min = 0.)).clip(max = 1.) * 255.
    # weight red channel
    c1 = 298.; c2 = 10.
    r_weight = np.arctan((ir - c1) / c2 * np.pi) / np.pi + 0.5
    red = r_weight * red + ((1. - r_weight) * (green > blue))
    # scale back to 0-1
    red /= red.max()
    green /= green.max()
    blue /= blue.max()

    # colour space transformation
    red2   = np.zeros(red.shape)
    green2 = np.zeros(red.shape)
    blue2  = np.zeros(red.shape)
    alpha  = np.ones(red.shape)
    # histogram for rescaling lightness - legacy from IDL code
    cumHist = np.array([0.00000, 0.141667, 0.245500, 0.311333, 0.381733, 0.446767, 0.507567, 
                        0.563467, 0.615867, 0.665867, 0.715867, 0.765067, 0.817533, 0.868433,
                        0.912533, 0.955433, 1.00000])
    print "transforming colour space"
    start_time = time.time()
    for i in range(red.shape[0]):
        for j in range(red.shape[1]):
            h, l, s = colorsys.rgb_to_hls(red[i, j], green[i, j], blue[i, j])
            h = (h + 24. / 360) % 1.
            l = np.interp(l*16, range(17), cumHist)
            red2[i, j], green2[i, j], blue2[i, j] = colorsys.hls_to_rgb(h, l, s**1.1)
            if np.isnan(red2[i, j]) or np.isnan(green2[i, j]) or np.isnan(blue2[i, j]):
                alpha[i, j] = 0. 
    elapsed_time = round(time.time() - start_time, 1)
    print "...done after " + str(elapsed_time) + " seconds"

    # further empirical scaling of RGB values
    # sea is blue, not red
    c1 = 0.8; c2 = 1. / c1; c3 = 1. - c2
    red2   = (((c3 + c2 * (red2  **c1)) * np.arctan(red2   * 100.) / np.pi * 2.)).clip(min = 0., max = 1.)
    c1 = 1.0; c2 = 1. / c1; c3 = 1. - c2
    blue2  = (((c3 + c2 * (blue2 **c1)) * np.arctan(blue2  * 100.) / np.pi * 2.)).clip(min = 0., max = 1.)
    c1 = 1.0; c2 = 1. / c1; c3 = 1. - c2
    green2 = (((c3 + c2 * (green2**c1)) * np.arctan(green2 * 100.) / np.pi * 2.)).clip(min = 0., max = 1.)

    # normalize to max = 1
    red2 /= np.nanmax(red2)
    green2 /= np.nanmax(green2)
    blue2 /= np.nanmax(blue2)

    # fix for matplotlib bug: exclude last row and column so that colourTuple agrees with coordinates
    red2   = red2  [0:-1, 0:-1]
    green2 = green2[0:-1, 0:-1]
    blue2  = blue2 [0:-1, 0:-1]
    alpha  = alpha [0:-1, 0:-1]

    # return colour array
    out = np.array([red2.flatten(), green2.flatten(), blue2.flatten(), alpha.flatten()]).transpose()
    return out

def plotRGB(figureName, colourTuple, lat, lon, dummy,
            lat_0, lon_0,
            width = 5000000, height = 5000000,
            res = 'l',
            llcrnrlat = 45, urcrnrlat = 80, llcrnrlon = -180, urcrnrlon = -150,
            cmin = None, cmax = None,
            mask = None):    

    # define map projection
    map = Basemap(width = width, height = height,
                  resolution='l', projection = 'stere',
                  lat_ts = lat_0, lat_0 = lat_0, lon_0 = lon_0)

    # define figure size
    fig1 = plt.figure(figsize = (10, 10))
    # define # subplots
    ax = fig1.add_subplot(111)

    # draw coasts and fill continents
    map.drawcoastlines(linewidth = 0.5)
    fillContinents = map.fillcontinents(color='#C0C0C0', lake_color='#7093DB', zorder = 0)
    
    # draw grid lines
    gridSpacing = 5.
    map.drawparallels(
        np.arange(llcrnrlat, urcrnrlat, gridSpacing),
        color = 'black', linewidth = 0.5,
        labels=[True, False, False, False])
    map.drawmeridians(
        np.arange(-180, 0, 10.),
        color = '0.25', linewidth = 0.5,
        labels=[False, False, False, True])
        
    # map.pcolormesh(lon, lat, x[:,:,0], latlon = True, linewidth = 0.05, clip_on = True)
    mesh = map.pcolormesh(lon, lat, dummy, color = colourTuple, latlon = True, linewidth = 0.05, clip_on = True)
    mesh.set_array(None)
    plt.savefig(figureName, bbox_inches='tight')

def plotRGBMulti(figureName, colourTuple,
            lat, lon, dummy,
            lat_0, lon_0,
            width = 5000000, height = 5000000,
            res = 'l',
            llcrnrlat = 45, urcrnrlat = 80, llcrnrlon = -180, urcrnrlon = -150,
            cmin = None, cmax = None,
            mask = None):    

    nPlots = colourTuple.shape

    # define map projection
    map = Basemap(width = width, height = height,
                  resolution='l', projection = 'stere',
                  lat_ts = lat_0, lat_0 = lat_0, lon_0 = lon_0)

    # define figure size
    fig1 = plt.figure(figsize = (20, 10))

    # FIGURE 1

    # define # subplots
    ax = fig1.add_subplot(121)

    # draw coasts and fill continents
    map.drawcoastlines(linewidth = 0.5)
    fillContinents = map.fillcontinents(color='#C0C0C0', lake_color='#7093DB', zorder = 0)
    
    # draw grid lines
    gridSpacing = 5.
    map.drawparallels(
        np.arange(llcrnrlat, urcrnrlat, gridSpacing),
        color = 'black', linewidth = 0.5,
        labels=[True, False, False, False])
    map.drawmeridians(
        np.arange(-180, 0, 10.),
        color = '0.25', linewidth = 0.5,
        labels=[False, False, False, True])
        
    # map.pcolormesh(lon, lat, x[:,:,0], latlon = True, linewidth = 0.05, clip_on = True)
    mesh = map.pcolormesh(lon, lat, dummy, color = colourTuple[:,:,0], latlon = True, linewidth = 0.05, clip_on = True)
    mesh.set_array(None)
    
    # FIGURE 2

    # define # subplots
    ax = fig1.add_subplot(122)

    # draw coasts and fill continents
    map.drawcoastlines(linewidth = 0.5)
    fillContinents = map.fillcontinents(color='#C0C0C0', lake_color='#7093DB', zorder = 0)
    
    # draw grid lines
    gridSpacing = 5.
    map.drawparallels(
        np.arange(llcrnrlat, urcrnrlat, gridSpacing),
        color = 'black', linewidth = 0.5,
        labels=[True, False, False, False])
    map.drawmeridians(
        np.arange(-180, 0, 10.),
        color = '0.25', linewidth = 0.5,
        labels=[False, False, False, True])
        
    # map.pcolormesh(lon, lat, x[:,:,0], latlon = True, linewidth = 0.05, clip_on = True)
    mesh = map.pcolormesh(lon, lat, dummy, color = colourTuple[:,:,1], latlon = True, linewidth = 0.05, clip_on = True)
    mesh.set_array(None)
    
    plt.savefig(figureName, bbox_inches='tight')
