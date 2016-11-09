from __future__ import division
from analyseCCI import CCI
from mpl_toolkits.basemap import Basemap, cm, latlon_default
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
from scipy import spatial
import sys
import numpy.ma as ma
from netCDF4 import Dataset
import colorsys
import time
from math import radians, cos, sin, asin, sqrt
from pandas import DataFrame, Index
import copy

def writeCCI(path, data, targetGrid, primary, platform = "N18"):
        
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

        phase = ncOut.createVariable("phase","f8",("along_track", "across_track"))
        phase[:,:] = data[:,:,13]

        cldtype = ncOut.createVariable("cldtype", "f8", ("along_track", "across_track"))
        cldtype[:, :] = data[:, :, 14]

    else:

        if platform == "N18":
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

        elif platform == "MYD":

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

        elif platform == "ENV":

            reflectance1 = ncOut.createVariable("reflectance_in_channel_no_1","f8",("along_track", "across_track"))
            reflectance1[:,:] = data[:,:,0]

            reflectance2 = ncOut.createVariable("reflectance_in_channel_no_2","f8",("along_track", "across_track"))
            reflectance2[:,:] = data[:,:,1]

            reflectance6 = ncOut.createVariable("reflectance_in_channel_no_3","f8",("along_track", "across_track"))
            reflectance6[:,:] = data[:,:,2]

            bt20 = ncOut.createVariable("brightness_temperature_in_channel_no_4","f8",("along_track", "across_track"))
            bt20[:,:] = data[:,:,3]

            bt31 = ncOut.createVariable("brightness_temperature_in_channel_no_5","f8",("along_track", "across_track"))
            bt31[:,:] = data[:,:,4]
        
            bt32 = ncOut.createVariable("brightness_temperature_in_channel_no_6","f8",("along_track", "across_track"))
            bt32[:,:] = data[:,:,5]
            
        else:
            print "ERROR in writeCCI: platform must be one of [N18, MYD, ENV]"
            sys.exit()

    ncOut.close()

def greatCircle(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return c * r

def resampleCCI(sourceCCI, targetGrid, sensor, lat_in = None, lon_in = None):
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
                     "phase", "cldtype"]
    else:
        if sensor == "N18":
            variables = ["albedo_in_channel_no_1", "albedo_in_channel_no_2", "albedo_in_channel_no_3",
                         "reflectance_in_channel_no_1", "reflectance_in_channel_no_2", "reflectance_in_channel_no_3",
                         "brightness_temperature_in_channel_no_4", "brightness_temperature_in_channel_no_5",
                         "brightness_temperature_in_channel_no_6"]    
        elif sensor == "MYD":
            variables = ["albedo_in_channel_no_1", "albedo_in_channel_no_2", "albedo_in_channel_no_6",
                         "reflectance_in_channel_no_1", "reflectance_in_channel_no_2", "reflectance_in_channel_no_6",
                         "brightness_temperature_in_channel_no_20", "brightness_temperature_in_channel_no_31",
                         "brightness_temperature_in_channel_no_32"]    
        elif sensor == "ENV":
            variables = ["reflec_nadir_0670", "reflec_nadir_0870", "reflec_nadir_1600",
                         "btemp_nadir_0370", "btemp_nadir_1100", "btemp_nadir_1200"]    
        else:
            print "ERROR in resampleCCI: sensor needs to be one of [N18, MYD, ENV]."

    if primary:
        sourceValues = np.zeros((np.prod(sourceCCI.time.shape), len(variables)))
    else:
        if sensor is "ENV":
            sourceValues = np.zeros((np.prod(sourceCCI.lat.shape), len(variables)))            
        else:
            sourceValues = np.zeros((np.prod(sourceCCI.cot_ap.shape), len(variables)))
    for var_i, var in enumerate(variables):
        sourceValues[:, var_i] = getattr(sourceCCI, var).ravel()
                 
    returnValues = np.zeros((targetPoints.shape[0], len(variables)))
    returnCount  = np.zeros((targetPoints.shape[0], len(variables)))

    i = 0
    for lon, lat in sourcePoints:
        nn = tree.query((lon, lat), k=1)
        targetIndex = nn[1]
        for var_i, var in enumerate(variables):
            """if data are not fill values"""
            if sourceValues[i, var_i] > -32767.0:
                """add value to box and increment counter"""
                """for phase, only average non-zero values"""
                if var is 'phase' and sourceValues[i, var_i] >= 1.:
                    returnCount[targetIndex, var_i] += 1
                    returnValues[targetIndex, var_i] += sourceValues[i, var_i]
                elif var is 'phase' and sourceValues[i, var_i] < 1.:
                    pass
                else:
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

def intToBin(x):
    return str(bin(x)[2:])

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
    elif platform is "ENV":
        ch4 = "brightness_temperature_in_channel_no_4"
        ch5 = "brightness_temperature_in_channel_no_5"        
    else:
        print "Error: select one of [MYD, N18] as input platform."
        sys.exit()
        
    # get variables if attributes not available
    if not hasattr(primaryData, "solar_zenith_view_no1"):
        primaryData.getAllVariables()
    if not hasattr(secondaryData, "reflectance_in_channel_no_1"):
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
        np.arange(-180, 180, 10.),
        color = '0.25', linewidth = 0.5,
        labels=[False, False, False, True])
        
    # map.pcolormesh(lon, lat, x[:,:,0], latlon = True, linewidth = 0.05, clip_on = True)
    mesh = map.pcolormesh(lon, lat, dummy, color = colourTuple, latlon = True, linewidth = 0.05, clip_on = True)
    mesh.set_array(None)
    plt.savefig(figureName, bbox_inches='tight')

def plotRGBMulti(figureName, colourTuple,
            lat, lon, dummy,
            latCalipso, lonCalipso,
            lat_0, lon_0,
            width = 5000000, height = 5000000,
            res = 'l',
            llcrnrlat = 45, urcrnrlat = 80, llcrnrlon = -180, urcrnrlon = -150,
            cmin = None, cmax = None,
            mask = None):    

    nPlots = colourTuple.shape[2]

    # define map projection
    map = Basemap(width = width, height = height,
                  resolution='l', projection = 'stere',
                  lat_ts = lat_0, lat_0 = lat_0, lon_0 = lon_0)

    # define calipso lat, lon
    x, y = np.reshape(lonCalipso, (4208, 1)), np.reshape(latCalipso, (4208, 1))
    xCalipso, yCalipso = map(x, y)

    # define figure size
    figWidth = 10 * nPlots
    fig1 = plt.figure(figsize = (figWidth, 10))

    # define # subplots
    nRows = 1
    nCols = nPlots
    
    lats = [58.5, 78, 83.5, 63]
    lons = [-121, -73, -63., -128.5]

    for i in range(nPlots):
        ax = fig1.add_subplot(nRows, nCols, i + 1)

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
            np.arange(-180, 180, 10.),
            color = '0.25', linewidth = 0.5,
            labels=[False, False, False, True])
        
        # map.pcolormesh(lon, lat, x[:,:,0], latlon = True, linewidth = 0.05, clip_on = True)
        mesh = map.pcolormesh(lon, lat, dummy, color = colourTuple[:,:,i], latlon = True, linewidth = 0.05, clip_on = True)
        mesh.set_array(None)
        mesh = map.plot(xCalipso, yCalipso, color = 'red', linestyle="dashed", linewidth=2.)
        draw_screen_poly(lats, lons, map)
    plt.savefig(figureName, bbox_inches='tight')
  
def mergeGranules(path1, path2, outPath):
    """Take two spatially adjacent MODIS granule files, read each variable, merge them, and write to new file"""
    print "Merging MODIS granules."
    
    # create CCI objects
    data1 = CCI(path1)
    data2 = CCI(path2)
    
    # get all variables
    data1.getAllVariables()
    data2.getAllVariables()
    
    # create output CCI object 
    data3 = CCI(path1)
    
    # create output NetCDF file
    ncOut = Dataset(outPath, "w", format="NETCDF4")
    
    # merge input variables and write to output variables
    i = 0
    for vName in iter(data1.dataset.variables):
        i += 1
        print vName
        v1 = getattr(data1, vName)
        v2 = getattr(data2, vName)
        
        vMerge = np.concatenate((v1,v2))
        setattr(data3, vName, vMerge)
        
        if i is 1:
            along_track = ncOut.createDimension("along_track", vMerge.shape[0])
            across_track = ncOut.createDimension("across_track", vMerge.shape[1])
    
        foo = ncOut.createVariable(vName,"f8",("along_track", "across_track"))
        foo[:,:] = getattr(data3, vName)        
    
    ncOut.close()

def draw_screen_poly(lats, lons, m):
    x, y = m(lons, lats)
    xy = zip(x,y)
    poly = Polygon(xy, color='red', alpha=1., fill=False, linewidth=2.)
    plt.gca().add_patch(poly)
    
def collocateCciAndCalipso(cci, calipso, maxDistance):
    print "collocating CCI and Calipso"
    cciLon = np.ma.compressed(cci.lon)
    cciLat = np.ma.compressed(cci.lat)
    cciCot = np.ma.compressed(cci.cot)
    cciCtt = np.ma.compressed(cci.ctt_corrected)
    cciCtp = np.ma.compressed(cci.ctp_corrected)
    cciCph = np.ma.compressed(cci.phase)
    cciCty = np.ma.compressed(cci.cldtype)

    calLon = calipso.get('lon')
    calLat = calipso.get('lat')
    calCod = calipso.get('cod')
    calCodLay = calipso.get('codLayered')
    calCtt = calipso.get('ctt')
    calCtp = calipso.get('ctp')
    calFcf = calipso.get('fcf')

    """ for each calipso pixel lat/lon,"""
    print "building source points"
    sourcePoints = np.asarray(zip(calLon, calLat))
    """ find closest CCI grid lat/lon"""
    print "building target points"
    targetPoints = np.asarray(zip(cciLon, cciLat))
    """ using a KDTree (which is built here)."""
    print "building KD trees"
    tree = spatial.cKDTree(targetPoints)

    """initialise collocated output variables"""
    colDist = []
    colCot = []
    colCtt = []
    colCtp = []
    colCph = []
    colCty = []
    colLatCalipso = []
    colLatCalipso1 = []
    colLatCalipso2 = []
    colLonCalipso = []
    colCodCalipso = []
    colCTP0Calipso = []
    colCTT0Calipso = []
    colCTP1Calipso = []
    colCTT1Calipso = []
    colCTP2Calipso = []
    colCTT2Calipso = []
    colPhase0Calipso = []
    colPhase1Calipso = []
    colPhase2Calipso = []
    colType0Calipso = []
    colType1Calipso = []
    colType2Calipso = []

    """ Loop over all Calipso lat/lons, """
    i = 0
    for lon, lat in sourcePoints:
        firstLayerFound = False
        secondLayerFound = False
        thirdLayerFound = False
        """ get the nearest CCI neighbour for each Calipso pixel, """
        nn = tree.query((lon, lat), k=1)
        """ extract its index in the flattened CCI lat/lons, """
        targetIndex = nn[1]
        lonNN = cciLon[targetIndex]
        latNN = cciLat[targetIndex]
        """ and calculate the great circle distance between Calipso lat/lon and nearest neighbour CCI lat/lon."""
        calipsoToBox = greatCircle(lon, lat, lonNN, latNN)
        """If this distance is smaller than the maximum possible distance, """
        if calipsoToBox < maxDistance:
            colCodCalipsoSum = 0.
            """" get Calipso feature flag, looping over atmosphere layers"""
            for l in range(calFcf.shape[1]):
                flagInt = int(calFcf[i, l])
                flagLength = flagInt.bit_length()
                flagBin = intToBin(calFcf[i, l])
                """ and extract cloud mask information at this layer."""
                colCmaskCalipso = int(flagBin[flagLength - 3:flagLength], 2)  # 0=invalid,1=clear,2=cloud
                """If Calipso says there is a cloud, """
                if colCmaskCalipso is 2:
                    """ sum up the layer optical depth until a threshold value is reached where we'll say the cloud top is."""
                    colCodCalipsoSum += calCodLay[i, l]
                    """If the Cod threshold has been exceeded, """
                    if colCodCalipsoSum > 0. and not firstLayerFound:
                        firstLayerFound = True
                        """add data to collocated variables"""
                        colDist.append(calipsoToBox)
                        colCot.append(cciCot[targetIndex])
                        colCtt.append(cciCtt[targetIndex])
                        colCtp.append(cciCtp[targetIndex])
                        colCph.append(cciCph[targetIndex])
                        colCty.append(cciCty[targetIndex])
                        colLatCalipso.append(lat)
                        colLonCalipso.append(lon)
                        colCodCalipso.append(calCod[i])
                        """get the phase and type"""
                        """PHASE: 1=ice,2=water,3=ice"""
                        colPhase0Calipso.append(int(flagBin[flagLength - 7:flagLength - 5], 2))
                        """TYPE: 0=low transp,1=low opaque,2=stratoc,3=low broken cum.,4=altocum.,5=altostr.,6=cirrus,7=deep conv."""
                        colType0Calipso.append(int(flagBin[flagLength - 12:flagLength - 9], 2))
                        colCTP0Calipso.append(calCtp[i, l])
                        colCTT0Calipso.append(calCtt[i, l] + 273.15)
                    if colCodCalipsoSum > 0.15 and not secondLayerFound:
                        """add data to collocated variables"""
                        """get the phase and type"""
                        secondLayerFound = True
                        colPhase1Calipso.append(int(flagBin[flagLength - 7:flagLength - 5], 2))  # 1=ice,2=water,3=ice
                        colType1Calipso.append(int(flagBin[flagLength - 12:flagLength - 9], 2))  # 0=low transp,1=low opaque,2=stratoc,3=low broken cum.,4=altocum.,5=altostr.,6=cirrus,7=deep conv.
                        colCTP1Calipso.append(calCtp[i, l])
                        colCTT1Calipso.append(calCtt[i, l] + 273.15)
                        colLatCalipso1.append(lat)
                    if colCodCalipsoSum > 1. and not thirdLayerFound:
                        """add data to collocated variables"""
                        """get the phase and type"""
                        thirdLayerFound = True
                        colPhase2Calipso.append(int(flagBin[flagLength - 7:flagLength - 5], 2))  # 1=ice,2=water,3=ice
                        colType2Calipso.append(int(flagBin[flagLength - 12:flagLength - 9], 2))  # 0=low transp,1=low opaque,2=stratoc,3=low broken cum.,4=altocum.,5=altostr.,6=cirrus,7=deep conv.
                        colCTP2Calipso.append(calCtp[i, l])
                        colCTT2Calipso.append(calCtt[i, l] + 273.15)
                        colLatCalipso2.append(lat)
                """if no additional cloud layers were found, add nan to layer variables"""
                if l == (calFcf.shape[1] - 1):
                    if not firstLayerFound:
                        colPhase0Calipso.append(np.nan)
                        colType0Calipso.append(np.nan)
                        colCTP0Calipso.append(np.nan)
                        colCTT0Calipso.append(np.nan)
                        colLatCalipso.append(lat)
                    if not secondLayerFound:
                        colPhase1Calipso.append(np.nan)
                        colType1Calipso.append(np.nan)
                        colCTP1Calipso.append(np.nan)
                        colCTT1Calipso.append(np.nan)
                        colLatCalipso1.append(lat)
                    if not thirdLayerFound:
                        colPhase2Calipso.append(np.nan)
                        colType2Calipso.append(np.nan)
                        colCTP2Calipso.append(np.nan)
                        colCTT2Calipso.append(np.nan)
                        colLatCalipso2.append(lat)

        i += 1

    """output variables should be numpy arrays and not dictionaries"""
    colCot = np.array(colCot)
    colCot = np.ma.masked_greater(colCot, 1000.)
    colCodCalipso = np.array(colCodCalipso)
    colCtt = np.array(colCtt)
    colCtt = np.ma.masked_greater(colCtt, 1000.)
    colCTTCalipso0 = np.array(colCTT0Calipso)
    colCTTCalipso1 = np.array(colCTT1Calipso)
    colCTTCalipso2 = np.array(colCTT2Calipso)
    colCtp = np.array(colCtp)
    colCtp = np.ma.masked_greater(colCtp, 10000.)
    colCTPCalipso0 = np.array(colCTP0Calipso)
    colCTPCalipso1 = np.array(colCTP1Calipso)
    colCTPCalipso2 = np.array(colCTP2Calipso)
    colLatCalipso = np.array(colLatCalipso)
    colCph = np.array(colCph)
    colCty = np.array(colCty)
    colPhase0Calipso = np.array(colPhase0Calipso)
    colPhase1Calipso = np.array(colPhase1Calipso)
    colPhase2Calipso = np.array(colPhase2Calipso)
    colType0Calipso = np.array(colType0Calipso)
    colType1Calipso = np.array(colType1Calipso)
    colType2Calipso = np.array(colType2Calipso)
    out = {'cciCot': colCot, 'cciCtt': colCtt, 'cciCtp': colCtp, 'cciCph': colCph, 'cciCty': colCty,
           'calipsoCtt0': colCTTCalipso0, 'calipsoCtp0': colCTPCalipso0,
           'calipsoPhase0': colPhase0Calipso, 'calipsoType0': colType0Calipso,
           'calipsoCtt1': colCTTCalipso1, 'calipsoCtp1': colCTPCalipso1,
           'calipsoPhase1': colPhase1Calipso, 'calipsoType1': colType1Calipso,
           'calipsoCtt2': colCTTCalipso2, 'calipsoCtp2': colCTPCalipso2,
           'calipsoPhase2': colPhase2Calipso, 'calipsoType2': colType2Calipso,
           'calipsoLat0': colLatCalipso, 'calipsoLat1': colLatCalipso1, 'calipsoLat2': colLatCalipso2,
           'calipsoCOD': colCodCalipso}
    return out

def plotCciCalipsoCollocation(collocateN18, collocateMYD, collocateENV, figuresDir):

    print "Plotting collocated data for Calipso and CCI."

    """First copy data dictionaries so that original values are preserved when manipulating data"""
    N18 = copy.deepcopy(collocateN18)
    MYD = copy.deepcopy(collocateMYD)
    ENV = copy.deepcopy(collocateENV)

    """The xaxis reference is Calipso's latitude"""
    plotLat0 = N18.get('calipsoLat0')
    plotLat1 = N18.get('calipsoLat1')
    plotLat2 = N18.get('calipsoLat2')
    """Figure settings."""
    fig = plt.figure(figsize=(15, 10))
    minX = plotLat0.min()
    maxX = plotLat0.max()
    """CTP"""
    ax = fig.add_subplot(2, 1, 1)  # 3 plots, vertically arranged
    ax.set_xlim([minX, maxX])
    plt.gca().invert_yaxis()
    ax.scatter(plotLat0, N18.get('cciCtp'), label="AVHRR")
    ax.scatter(plotLat0, MYD.get('cciCtp'), label="MODIS AQUA", c="g")
    ax.scatter(plotLat0, ENV.get('cciCtp'), label="AATSR", c="white")
    ax.scatter(plotLat2, N18.get('calipsoCtp2'), label="Calipso [COT > 1]", c="pink")
    ax.scatter(plotLat0, N18.get('calipsoCtp0'), label="Calipso [COT > 0]", c="r")
    ax.set_ylabel("CTP [hPa]")
    handles, labels = ax.get_legend_handles_labels()
    order=[0,1,2,4,3]
    labels = [ labels[i] for i in order]
    handles = [handles[i] for i in order]
    leg = ax.legend(handles=handles, labels=labels, loc=3, frameon=True, fancybox=True, fontsize=11)
    leg.get_frame().set_alpha(0.5)
    # """CTT"""
    # ax = fig.add_subplot(3, 1, 2)
    # plt.gca().invert_yaxis()
    # ax.set_xlim([minX, maxX])
    # ax.scatter(plotLat0, N18.get('cciCtt'))
    # ax.scatter(plotLat0, MYD.get('cciCtt'), c="g")
    # ax.scatter(plotLat0, ENV.get('cciCtt'), c="y")
    # ax.scatter(plotLat0, N18.get('calipsoCtt0'), c="r")
    # ax.scatter(plotLat1, N18.get('calipsoCtt1'), c="pink")
    # ax.set_ylabel("CTT [K]")
    """COT"""
    ax = fig.add_subplot(2, 1, 2)
    ax.set_ylim([0, 50])
    plt.gca().invert_yaxis()
    ax.set_xlim([minX, maxX])
    ax.scatter(plotLat0, N18.get('cciCot'))
    ax.scatter(plotLat0, MYD.get('cciCot'), c="g")
    ax.scatter(plotLat0, ENV.get('cciCot'), c="white")
    ax.scatter(plotLat0, N18.get('calipsoCOD'), c="r")
    ax.set_ylabel("COT")
    plt.xlabel("Latitude")
    """CLOUD PHASE TABLE"""
    idx = Index(range(1, 6))
    cphCal0 = np.array(N18.get('calipsoPhase0'))
    cphCal0 = cphCal0.astype(float)
    cphCal0[cphCal0 == 1] = 999.
    cphCal0[cphCal0 == 2] = 1.
    cphCal0[cphCal0 == 999.] = 2.
    cphCal0[cphCal0 == 3] = 2.
    cphCal0[cphCal0 == 0.] = 999. #np.nan
    # cphCal0[cphCal0 > 3.] = np.nan
    cphCal1 = np.array(N18.get('calipsoPhase1'))
    cphCal1 = cphCal1.astype(float)
    cphCal1[cphCal1 == 1] = 999.
    cphCal1[cphCal1 == 2] = 1.
    cphCal1[cphCal1 == 999.] = 2.
    cphCal1[cphCal1 == 3] = 2.
    cphCal1[cphCal1 == 0.] = 999. #np.nan
    # cphCal1[cphCal1 > 3.] = np.nan
    cphCal2 = np.array(N18.get('calipsoPhase2'))
    cphCal2 = cphCal2.astype(float)
    cphCal2[cphCal2 == 1] = 999.
    cphCal2[cphCal2 == 2] = 1.
    cphCal2[cphCal2 == 999.] = 2.
    cphCal2[cphCal2 == 3] = 2.
    cphCal2[cphCal2 == 0.] = 999. #np.nan
    # cphCal2[cphCal2 > 3.] = np.nan
    cty0 = np.array(N18.get('calipsoType0'))
    cty1 = np.array(N18.get('calipsoType1'))
    cty2 = np.array(N18.get('calipsoType2'))
    cphN18 = N18.get('cciCph')
    cphN18[cphN18 == 0.] = np.nan
    cphN18[cphN18 > 2.] = np.nan
    ctyN18 = np.round(N18.get('cciCty'), 0)
    ctyN18[ctyN18 == 0.] = np.nan
    cphMYD = MYD.get('cciCph')
    cphMYD[cphMYD == 0.] = np.nan
    cphMYD[cphMYD > 2.] = np.nan
    ctyMYD = str(np.round(MYD.get('cciCty'), 0))
    ctyMYD[ctyMYD == 0.] = '' #np.nan
    cphENV = ENV.get('cciCph')
    cphENV[cphENV == 0.] = np.nan
    cphENV[cphENV > 2.] = np.nan
    ctyENV = np.round(ENV.get('cciCty'), 0)
    ctyENV[ctyENV == 0.] = np.nan
    df = DataFrame([cphCal0, cphCal1, cphCal2, cphN18, cphMYD, cphENV])
    vals = np.around(df.values, 2)
    normal = plt.Normalize(np.nanmin(vals) - 1, np.nanmax(vals) + 1)
    cell_colours = plt.cm.RdBu(normal(vals))
    """cell colours scale from blue (phase = 1) to red (phase = 2)
        i.e. green is 0, and only red and blue depend on phase:
        red  = phase - 1.
        blue = 1. - red"""
    cell_colours[:,:, 0] = vals - 1.  # red
    cell_colours[:, :, 1] = 0.        # green
    cell_colours[:, :, 2] = 2. - vals # blue
    cell_colours[np.isnan(vals), 0:3] = 1.
    cell_colours[abs(cell_colours) > 900.] = 0.5
    cell_colours[cell_colours[0:3,:,0]==0.5, 1] = 0.5
    #foo = np.ones((1,120,4))
    #cell_colours = np.vstack((foo, cell_colours))
    row_labels = ["Calipso [COT > 0]", "Calipso [COT > 0.15]", "Calipso [COT > 1]", "AVHRR", "MODIS AQUA", "AATSR"]
    cell_text = np.chararray((len(row_labels), len(plotLat0)))
    cell_text[0,] = cty0
    cell_text[1,] = cty1
    cell_text[2,] = cty2
    cell_text[3,] = ctyN18
    cell_text[4,] = ctyMYD
    cell_text[5,] = ctyENV
    table = plt.table(cellText=cell_text,
                      rowLabels=row_labels,
                      cellColours=cell_colours,
                      bbox=[0., -0.45, 1., 0.3],
                      loc='bottom')
    table.set_fontsize(10)
    # iterate through cells of a table
    table_props = table.properties()
    table_cells = table_props['child_artists']
    for cell in table_cells:
        cell._text.set_color('white')
    table._cells[(1, 0)]._text.set_color('black')
    table.set_fontsize(10)
    """save figure"""
    plt.savefig(figuresDir + 'calipsoVsCci.png', bbox_inches='tight')

    
    

