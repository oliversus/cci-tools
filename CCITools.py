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
import os
from operator import itemgetter
import globals
import glob, re

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

        cot_unc = ncOut.createVariable("cot_uncertainty","f8",("along_track", "across_track"))
        cot_unc[:,:] = data[:,:,4]

        cer = ncOut.createVariable("cer","f8",("along_track", "across_track"))
        cer[:,:] = data[:,:,5]
        
        cwp = ncOut.createVariable("cwp","f8",("along_track", "across_track"))
        cwp[:,:] = data[:,:,6]
        
        ctp = ncOut.createVariable("ctp","f8",("along_track", "across_track"))
        ctp[:,:] = data[:,:,7]

        ctp_unc = ncOut.createVariable("ctp_uncertainty","f8",("along_track", "across_track"))
        ctp_unc[:,:] = data[:,:,8]
        
        ctp_corrected = ncOut.createVariable("ctp_corrected","f8",("along_track", "across_track"))
        ctp_corrected[:,:] = data[:,:,9]
        
        cth = ncOut.createVariable("cth","f8",("along_track", "across_track"))
        cth[:,:] = data[:,:,10]

        cth_unc = ncOut.createVariable("cth_uncertainty","f8",("along_track", "across_track"))
        cth_unc[:,:] = data[:,:,11]

        cth_corrected = ncOut.createVariable("cth_corrected","f8",("along_track", "across_track"))
        cth_corrected[:,:] = data[:,:,12]
        
        ctt = ncOut.createVariable("ctt","f8",("along_track", "across_track"))
        ctt[:,:] = data[:,:,13]

        ctt_uncertainty = ncOut.createVariable("ctt_uncertainty","f8",("along_track", "across_track"))
        ctt_uncertainty[:,:] = data[:,:,14]

        ctt_corrected = ncOut.createVariable("ctt_corrected","f8",("along_track", "across_track"))
        ctt_corrected[:,:] = data[:,:,15]

        cc_total = ncOut.createVariable("cc_total","f8",("along_track", "across_track"))
        cc_total[:,:] = data[:,:,16]

        phase = ncOut.createVariable("phase","f8",("along_track", "across_track"))
        phase[:,:] = data[:,:,17]

        cldtype = ncOut.createVariable("cldtype", "f8", ("along_track", "across_track"))
        cldtype[:, :] = data[:, :, 18]

        nisemask = ncOut.createVariable("nisemask", "f8", ("along_track", "across_track"))
        nisemask[:, :] = data[:, :, 19]

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

def resampleCCI(sourceCCI, targetGrid, sensor, maxDistance, lat_in = None, lon_in = None):
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
        variables = ["time", "solar_zenith_view_no1", "satellite_zenith_view_no1", "cot", "cot_uncertainty", "cer", "cwp",
                     "ctp", "ctp_uncertainty", "ctp_corrected", "cth", "cth_uncertainty", "cth_corrected", "ctt", "ctt_uncertainty", "ctt_corrected", "cc_total",
                     "phase", "cldtype", "nisemask"]
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
    returnCldType = [[] for i in range(targetPoints.shape[0])]
    targetLon = targetGrid.lon.flatten()
    targetLat = targetGrid.lat.flatten()

    i = 0
    for lon, lat in sourcePoints:
        nn = tree.query((lon, lat), k=1)
        targetIndex = nn[1]
        """Check with target grid lat/lons, """
        lonNN = targetLon[targetIndex]
        latNN = targetLat[targetIndex]
        """ that the great circle distance between CCI L2 lat/lon and target grid box lat/lon is lower than a threshold value (i.e. L2 pixel is within grid box)."""
        L2ToBox = greatCircle(lon, lat, lonNN, latNN)
        if L2ToBox < maxDistance:
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
                    elif var is 'cldtype' and sourceValues[i, var_i] >= 1.:
                        returnCldType[targetIndex].append(sourceValues[i, var_i])
                    elif var is 'cldtype' and sourceValues[i, var_i] < 1.:
                        pass
                    else:
                        returnCount [targetIndex, var_i] += 1
                        returnValues[targetIndex, var_i] += sourceValues[i, var_i]
        i += 1
        if i % 250000 == 0:
            print ('{0}\r'.format(round(i / sourcePoints.shape[0] * 100, 1)))   

    returnCountMasked = ma.masked_values(returnCount, 0.)
    out = returnValues / returnCountMasked

    if primary:
        for i in range(len(returnCldType)):
            phase = round(out[i, variables.index('phase')], 0)
            cldType = np.array(returnCldType[i]).astype(int)
            if phase == 1.:
                cldType = cldType[cldType < 5]
                out[i, variables.index('cldtype')] = randomMode(cldType)
            elif phase == 2:
                cldType = cldType[cldType >= 5]
                out[i, variables.index('cldtype')] = randomMode(cldType)
            else:
                pass

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

    diff = False
    # output figure name:
    if platform == "NOAA18":
        CCI_filename = os.path.basename(priN18.getPath())
    elif platform == "MYD":
        CCI_filename = os.path.basename(priMYD.getPath())
    elif platform == "ENV":
        CCI_filename = os.path.basename(priENV.getPath())
    else:
        CCI_filename = os.path.basename(priN18.getPath())
        diff = True

    sep = "_"; suffix = ".png"
    foo = list(itemgetter(3, 4, 6, 7)(CCI_filename.replace(".nc", "_nc").split("_")).__add__((variable,)))
    if diff:
        foo[0] = platform
    figure_name = globals.figuresDir + sep.join(foo) + suffix

    # define figure size
    #fig1 = plt.figure(figsize=(10, 10))

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
                 cmin = cmin, cmax = cmax,
                 mask = mask)
    plt.savefig(figure_name, bbox_inches='tight')

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
    fillContinents = map.fillcontinents(color='#C0C0C0', lake_color='#7093DB', zorder=0)
    
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
    mesh = map.pcolormesh(lon, lat, dummy, color=colourTuple, latlon=True, linewidth=0.05, clip_on=True)
    mesh.set_array(None)
    plt.savefig(figureName, bbox_inches='tight')

def plotRGBMulti(figureName, colourTuple,
            lat, lon, poly_lats, poly_lons, dummy,
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
    x, y = np.reshape(lonCalipso, (lonCalipso.shape[0], 1)), np.reshape(latCalipso, (latCalipso.shape[0], 1))
    xCalipso, yCalipso = map(x, y)

    # define figure size
    figWidth = 10 * nPlots
    fig1 = plt.figure(figsize = (figWidth, 10))

    # define # subplots
    nRows = 1
    nCols = nPlots

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
        draw_screen_poly(poly_lats, poly_lons, map)
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
    
def collocateCciAndCalipso(cci, calipso, maxDistance, corrected):
    print "collocating CCI and Calipso"
    cciLon = np.ma.compressed(cci.lon)
    cciLat = np.ma.compressed(cci.lat)
    cciCot = np.ma.compressed(cci.cot)
    cciCotUnc = np.ma.compressed(cci.cot_uncertainty)
    if corrected:
        cciCtt = np.ma.compressed(cci.ctt_corrected)
        cciCtp = np.ma.compressed(cci.ctp_corrected)
    else:
        cciCtt = np.ma.compressed(cci.ctt)
        cciCtp = np.ma.compressed(cci.ctp)
    cciCtpUnc = np.ma.compressed(cci.ctp_uncertainty)
    cciCph = np.ma.compressed(cci.phase)
    cciCty = np.ma.compressed(cci.cldtype)
    cciNise = np.ma.compressed(cci.nisemask)

    calLon = calipso.get('lon')
    calLat = calipso.get('lat')
    calCod = calipso.get('cod')
    calCodLay = calipso.get('codLayered')
    calCtt = calipso.get('ctt')
    calCth = calipso.get('cth')
    calCtp = calipso.get('ctp')
    calCtpBot = calipso.get('ctpBot')
    calFcf = calipso.get('fcf')
    calIce = calipso.get('ice')
    calTop = calipso.get('top')
    calTyp = calipso.get('typ')

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
    colCotUnc = []
    colCtt = []
    colCtp = []
    colCtpUnc = []
    colCph = []
    colCty = []
    colNise = []
    colLatCalipso = []
    colLatCalipso1 = []
    colLatCalipso2 = []
    colLonCalipso = []
    colCodCalipso = []
    colCTP0Calipso = []
    colCTPBot0Calipso = []
    colCTT0Calipso = []
    colCTP1Calipso = []
    colCTPBot1Calipso = []
    colCTT1Calipso = []
    colCTP2Calipso = []
    colCTPBot2Calipso = []
    colCTT2Calipso = []
    colCTHCalipso = []
    colPhase0Calipso = []
    colPhase1Calipso = []
    colPhase2Calipso = []
    colType0Calipso = []
    colType1Calipso = []
    colType2Calipso = []
    colIceCalipso = []
    colTopCalipso = []
    colTypCalipso = []

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
            colIceCalipso.append(calIce[i, 0])
            colTopCalipso.append(calTop[i, 0])
            colTypCalipso.append(calTyp[i, 0])
            colCodCalipsoSum = 0.
            """" get Calipso feature flag, looping over atmosphere layers"""
            colDist.append(calipsoToBox)
            colCot.append(cciCot[targetIndex])
            colCotUnc.append(cciCotUnc[targetIndex])
            colCtt.append(cciCtt[targetIndex])
            colCtp.append(cciCtp[targetIndex])
            colCtpUnc.append(cciCtpUnc[targetIndex])
            colCph.append(cciCph[targetIndex])
            colCty.append(cciCty[targetIndex])
            colNise.append(cciNise[targetIndex])
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
                        colCodCalipso.append(calCod[i])
                        """get the phase and type"""
                        """PHASE: 1=ice,2=water,3=ice"""
                        colPhase0Calipso.append(int(flagBin[flagLength - 7:flagLength - 5], 2))
                        """TYPE: 0=low transp,1=low opaque,2=stratoc,3=low broken cum.,4=altocum.,5=altostr.,6=cirrus,7=deep conv."""
                        colType0Calipso.append(int(flagBin[flagLength - 12:flagLength - 9], 2))
                        colCTP0Calipso.append(calCtp[i, l])
                        colCTPBot0Calipso.append(calCtpBot[i, l])
                        colCTT0Calipso.append(calCtt[i, l] + 273.15)
                        colCTHCalipso.append(calCth[i, l])
                    if colCodCalipsoSum > 0.15 and not secondLayerFound:
                        """add data to collocated variables"""
                        """get the phase and type"""
                        secondLayerFound = True
                        colPhase1Calipso.append(int(flagBin[flagLength - 7:flagLength - 5], 2))  # 1=ice,2=water,3=ice
                        colType1Calipso.append(int(flagBin[flagLength - 12:flagLength - 9], 2))  # 0=low transp,1=low opaque,2=stratoc,3=low broken cum.,4=altocum.,5=altostr.,6=cirrus,7=deep conv.
                        colCTP1Calipso.append(calCtp[i, l])
                        colCTPBot1Calipso.append(calCtpBot[i, l])
                        colCTT1Calipso.append(calCtt[i, l] + 273.15)
                        colLatCalipso1.append(lat)
                    if colCodCalipsoSum > 1. and not thirdLayerFound:
                        """add data to collocated variables"""
                        """get the phase and type"""
                        thirdLayerFound = True
                        colPhase2Calipso.append(int(flagBin[flagLength - 7:flagLength - 5], 2))  # 1=ice,2=water,3=ice
                        colType2Calipso.append(int(flagBin[flagLength - 12:flagLength - 9], 2))  # 0=low transp,1=low opaque,2=stratoc,3=low broken cum.,4=altocum.,5=altostr.,6=cirrus,7=deep conv.
                        colCTP2Calipso.append(calCtp[i, l])
                        colCTPBot2Calipso.append(calCtpBot[i, l])
                        colCTT2Calipso.append(calCtt[i, l] + 273.15)
                        colLatCalipso2.append(lat)

                if l == (calFcf.shape[1] - 1):
                    """if in last layer, get calipso lat/lon"""
                    colLatCalipso.append(lat)
                    colLonCalipso.append(lon)
                    """if no additional cloud layers were found, add nan to layer variables"""
                    if not firstLayerFound:
                        colPhase0Calipso.append(np.nan)
                        colType0Calipso.append(np.nan)
                        colCTP0Calipso.append(np.nan)
                        colCTPBot0Calipso.append(np.nan)
                        colCTT0Calipso.append(np.nan)
                        colCodCalipso.append(calCod[i])
                        colCTHCalipso.append(np.nan)
                    if not secondLayerFound:
                        colPhase1Calipso.append(np.nan)
                        colType1Calipso.append(np.nan)
                        colCTP1Calipso.append(np.nan)
                        colCTPBot1Calipso.append(np.nan)
                        colCTT1Calipso.append(np.nan)
                        colLatCalipso1.append(lat)
                    if not thirdLayerFound:
                        colPhase2Calipso.append(np.nan)
                        colType2Calipso.append(np.nan)
                        colCTP2Calipso.append(np.nan)
                        colCTPBot2Calipso.append(np.nan)
                        colCTT2Calipso.append(np.nan)
                        colLatCalipso2.append(lat)
        i += 1

    """output variables should be numpy arrays and not dictionaries"""
    colCot = np.array(colCot)
    colCot = np.ma.masked_greater(colCot, 1000.)
    colCotUnc = np.array(colCotUnc)
    colCotUnc = np.ma.masked_greater(colCotUnc, 1000.)
    colCodCalipso = np.array(colCodCalipso)
    colCtt = np.array(colCtt)
    colCtt = np.ma.masked_greater(colCtt, 1000.)
    colCTTCalipso0 = np.array(colCTT0Calipso)
    colCTTCalipso1 = np.array(colCTT1Calipso)
    colCTTCalipso2 = np.array(colCTT2Calipso)
    colCTHCalipso = np.array(colCTHCalipso)
    colCtp = np.array(colCtp)
    colCtp = np.ma.masked_greater(colCtp, 10000.)
    colCtpUnc = np.array(colCtpUnc)
    colCtpUnc = np.ma.masked_greater(colCtpUnc, 10000.)
    colCTPCalipso0 = np.array(colCTP0Calipso)
    colCTPCalipso1 = np.array(colCTP1Calipso)
    colCTPCalipso2 = np.array(colCTP2Calipso)
    colCTPBotCalipso0 = np.array(colCTPBot0Calipso)
    colCTPBotCalipso1 = np.array(colCTPBot1Calipso)
    colCTPBotCalipso2 = np.array(colCTPBot2Calipso)
    colLatCalipso = np.array(colLatCalipso)
    colCph = np.array(colCph)
    colCty = np.array(colCty)
    colNise = np.array(colNise)
    colPhase0Calipso = np.array(colPhase0Calipso)
    colPhase1Calipso = np.array(colPhase1Calipso)
    colPhase2Calipso = np.array(colPhase2Calipso)
    colType0Calipso = np.array(colType0Calipso)
    colType1Calipso = np.array(colType1Calipso)
    colType2Calipso = np.array(colType2Calipso)
    colIceCalipso = np.array(colIceCalipso)
    colTopCalipso = np.array(colTopCalipso)
    colTypCalipso = np.array(colTypCalipso)
    out = {'cciCot': colCot, 'cciCotUnc': colCotUnc, 'cciCtt': colCtt, 'cciCtp': colCtp, 'cciCtpUnc': colCtpUnc,
           'cciCph': colCph, 'cciCty': colCty, 'cciNise': colNise,
           'calipsoCtt0': colCTTCalipso0, 'calipsoCtp0': colCTPCalipso0, 'calipsoCtpBot0': colCTPBotCalipso0,
           'calipsoPhase0': colPhase0Calipso, 'calipsoType0': colType0Calipso,
           'calipsoCtt1': colCTTCalipso1, 'calipsoCtp1': colCTPCalipso1, 'calipsoCtpBot1': colCTPBotCalipso1,
           'calipsoPhase1': colPhase1Calipso, 'calipsoType1': colType1Calipso,
           'calipsoCtt2': colCTTCalipso2, 'calipsoCtp2': colCTPCalipso2, 'calipsoCtpBot2': colCTPBotCalipso2,
           'calipsoPhase2': colPhase2Calipso, 'calipsoType2': colType2Calipso, 'calipsoCth': colCTHCalipso,
           'calipsoLat0': colLatCalipso, 'calipsoLat1': colLatCalipso1, 'calipsoLat2': colLatCalipso2,
           'calipsoCOD': colCodCalipso, 'calipsoIce': colIceCalipso, 'calipsoTop': colTopCalipso, 'calipsoTyp': colTypCalipso}
    return out

def plotCciCalipsoCollocation(collocateN18, collocateMYD, collocateENV, figurePath, sceneTime, plotCot):

    print "Plotting collocated data for Calipso and CCI."

    """First copy data dictionaries so that original values are preserved when manipulating data"""
    N18 = copy.deepcopy(collocateN18)
    MYD = copy.deepcopy(collocateMYD)
    ENV = copy.deepcopy(collocateENV)

    """calculate CCI COT minus Calipso COD = deltaCOT"""
    """next step:
        plot only for those pixels where both say ice"""
    calipso_variable = 'calipsoCtp0'
    cc4cl_variable = 'cciCtp'
    fig = plt.figure(figsize=(20, 10))
    """loop over both phases"""
    # write_NA1 = True
    # write_SIB = True
    for i in range(1, 3):
        for j in range(0, 3):
            N18_phase = copy.deepcopy(collocateN18)
            MYD_phase = copy.deepcopy(collocateMYD)
            ENV_phase = copy.deepcopy(collocateENV)
            cphCal0 = correct_calipso_phase(N18_phase.get('calipsoPhase' + str(j)))
            cphN18 = np.round(correct_cci_phase(N18_phase.get('cciCph')))
            cphMYD = np.round(correct_cci_phase(MYD_phase.get('cciCph')))
            cphENV = np.round(correct_cci_phase(ENV_phase.get('cciCph')))
            """phase to be matched"""
            phase_to_match = i
            """is used to remove non-phase-matching pixels"""
            is_no_phase_match_N18 = np.logical_or(cphCal0 != phase_to_match, cphN18 != phase_to_match)
            is_no_phase_match_MYD = np.logical_or(cphCal0 != phase_to_match, cphMYD != phase_to_match)
            is_no_phase_match_ENV = np.logical_or(cphCal0 != phase_to_match, cphENV != phase_to_match)
            """then, calculate difference between CCI and calipso"""
            deltaCtpN18 = N18_phase.get(cc4cl_variable) - N18_phase.get(calipso_variable)
            deltaCtpMYD = MYD_phase.get(cc4cl_variable) - MYD_phase.get(calipso_variable)
            deltaCtpENV = ENV_phase.get(cc4cl_variable) - ENV_phase.get(calipso_variable)
            """before accounting for phase matches, calculate biases for paper statistics"""
            # if sceneTime is globals.NA1 and write_NA1:
            if sceneTime is globals.NA1:
                NA1_ctp_bias_part = round(np.nanmean(deltaCtpN18[N18_phase.get('calipsoLat0') > 73.7]), 2)
                globals.latex_variables['NA1_ctp_bias_part'] = NA1_ctp_bias_part
                #write = "NA1_ctp_bias_part=" + str(NA1_ctp_bias_part) + "\n"
                #latex_file.write(write)
                # write_NA1 = False
            if sceneTime is globals.SIB:
                N18_bias = np.nanmean(deltaCtpN18[N18_phase.get('calipsoLat0') > 74.])
                MYD_bias = np.nanmean(deltaCtpMYD[MYD_phase.get('calipsoLat0') > 74.])
                ENV_bias = np.nanmean(deltaCtpENV[ENV_phase.get('calipsoLat0') > 74.])
                SIB_ctp_bias_part = round(N18_bias, MYD_bias, ENV_bias, 2)
                globals.latex_variables['SIB_ctp_bias_part'] = SIB_ctp_bias_part
                # write = "SIB_ctp_bias_part=" + str(SIB_ctp_bias_part) + "\n"
                # latex_file.write(write)
                # write_SIB = False

            """now, account for phase matching"""
            deltaCtpN18[is_no_phase_match_N18] = np.nan
            deltaCtpMYD[is_no_phase_match_MYD] = np.nan
            deltaCtpENV[is_no_phase_match_ENV] = np.nan

            """and get retrieval uncertainties"""
            cotUncN18 = N18_phase.get(cc4cl_variable + 'Unc')
            cotUncN18[is_no_phase_match_N18] = np.nan
            cotUncMYD = MYD_phase.get(cc4cl_variable + 'Unc')
            cotUncMYD[is_no_phase_match_MYD] = np.nan
            cotUncENV = ENV_phase.get(cc4cl_variable + 'Unc')
            cotUncENV[is_no_phase_match_ENV] = np.nan
            """plot these values"""
            panel = 3 * i + j - 2
            ax = fig.add_subplot(2, 3, panel)
            ax.scatter(cotUncN18, deltaCtpN18, label="AVHRR")
            ax.scatter(cotUncMYD, deltaCtpMYD, label="MODIS AQUA", c="red")
            ax.scatter(cotUncENV, deltaCtpENV, label="AATSR", c="green")
            if phase_to_match == 1:
                phase = "water"
            else:
                phase = "ice"
            ax.set_xlabel("CC4CL CTP uncertainty (" + phase + ", layer " + str(j) + ")")
            ax.set_ylabel("CC4CL CTP - Calipso CTP")
            if panel is 1:
                leg = ax.legend(loc=2, frameon=True, fancybox=True, fontsize=11)
                leg.get_frame().set_alpha(0.5)
    figurePathUnc = os.path.dirname(figurePath) + "/" + "_".join(figurePath[figurePath.find("calipsoVsCci"):-1].split("_")[0:2]) + "_uncertainty.png"
    plt.savefig(figurePathUnc, bbox_inches='tight')


    """The xaxis reference is Calipso's latitude"""
    plotLat0 = N18.get('calipsoLat0')
    plotLat2 = N18.get('calipsoLat2')
    """Figure settings."""
    fig = plt.figure(figsize=(15, 10))
    minX = plotLat0.min()
    maxX = plotLat0.max()
    """CTP"""
    ax = fig.add_subplot(2, 1, 1)  # several plots, vertically arranged
    width = (max(plotLat2)-min(plotLat2)) / len(plotLat2) # bar width
    ax.set_xlim([minX - width * 0.5, maxX + width * 0.5])
    ax.set_ylim([100, 1000])
    plt.gca().invert_yaxis()
    if not plotCot:
        plot_topography(N18, plotLat0, ax)
        ax.set_xlabel("Latitude")
    alpha = 1
    ax.bar(plotLat2 - width * 0.5, N18.get('calipsoCtp0') - N18.get('calipsoCtpBot0'), bottom=N18.get('calipsoCtpBot0'), width=width,
           color="lavenderblush", alpha=alpha, label="Calipso [COT > 0]", zorder=3)
    ax.bar(plotLat2 - width * 0.5, N18.get('calipsoCtp1') - N18.get('calipsoCtpBot1'), bottom=N18.get('calipsoCtpBot1'), width=width,
           color="lightpink", alpha=alpha, label="Calipso [COT > 0.15]", zorder=3)
    ax.bar(plotLat2 - width * 0.5, N18.get('calipsoCtp2') - N18.get('calipsoCtpBot2'), bottom=N18.get('calipsoCtpBot2'), width=width,
           color="thistle", alpha=alpha, label="Calipso [COT > 1]", zorder=3)
    ax.scatter(plotLat0, ENV.get('cciCtp'), label="AATSR", c="orange", zorder=10)
    ax.scatter(plotLat0, MYD.get('cciCtp'), label="MODIS AQUA", c="aqua", zorder=10)
    ax.scatter(plotLat0, N18.get('cciCtp'), label="AVHRR", c="r", zorder=10)
    ax.set_ylabel("CTP [hPa]")
    handles, labels = ax.get_legend_handles_labels()
    order=[2,1,0,3,4,5]
    labels = [ labels[i] for i in order]
    handles = [handles[i] for i in order]
    if sceneTime == '07222058':
        loc = 2
    elif sceneTime == '07270810':
        loc = 1
    elif sceneTime == '07221915':
        loc = 1
    else:
        loc = 3
    leg = ax.legend(handles=handles, labels=labels, loc=loc, frameon=True, fancybox=True, fontsize=11)
    leg.get_frame().set_alpha(0.5)

    """COT"""
    if plotCot:
        ax = fig.add_subplot(2, 1, 2)
        ax.set_ylim([0, 50])
        plt.gca().invert_yaxis()
        ax.set_xlim([minX - width * 0.5, maxX + width * 0.5])
        plot_topography(N18, plotLat0, ax)
        ax.scatter(plotLat0, MYD.get('cciCot'), c="orange", zorder=10)
        ax.scatter(plotLat0, ENV.get('cciCot'), c="aqua", zorder=10)
        ax.scatter(plotLat0, N18.get('cciCot'), c="r", zorder=10)
        foo = np.where(N18.get('calipsoCOD')==0, np.nan, N18.get('calipsoCOD'))
        ax.scatter(plotLat0, foo, c="thistle", zorder=10)
        ax.set_ylabel("COT")
        ax.set_xlabel("Latitude")

    """CLOUD PHASE TABLE"""
    cphCal0 = correct_calipso_phase(N18.get('calipsoPhase0'))
    cphCal1 = correct_calipso_phase(N18.get('calipsoPhase1'))
    cphCal2 = correct_calipso_phase(N18.get('calipsoPhase2'))
    cty0 = np.array(N18.get('calipsoType0'))
    cty1 = np.array(N18.get('calipsoType1'))
    cty2 = np.array(N18.get('calipsoType2'))
    cphN18 = correct_cci_phase(N18.get('cciCph'))
    ctyN18 = np.round(N18.get('cciCty'), 0)
    ctyN18[ctyN18 == 0.] = np.nan
    cphMYD = correct_cci_phase(MYD.get('cciCph'))
    ctyMYD = np.round(MYD.get('cciCty'), 0)
    ctyMYD[ctyMYD == 0.] = np.nan
    cphENV = correct_cci_phase(ENV.get('cciCph'))
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
    cell_colours[:,:,3] = 0.8
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
    # iterate through cells of a table
    table_props = table.properties()
    table_cells = table_props['child_artists']
    for cell in table_cells:
        cell._text.set_color('white')
    table._cells[(0, -1)]._text.set_color('black')
    table._cells[(1, -1)]._text.set_color('black')
    table._cells[(2, -1)]._text.set_color('black')
    table._cells[(3, -1)]._text.set_color('black')
    table._cells[(4, -1)]._text.set_color('black')
    table._cells[(5, -1)]._text.set_color('black')
    """save figure"""
    plt.savefig(figurePath, bbox_inches='tight')

def randomMode(data, max=11):
    """After counting the number of occurrences of each integer within data"""
    counts = np.bincount(data, minlength=max)
    """check whether several maxima exist"""
    max_indexes = np.where(counts == counts.max())[0]
    if (len(max_indexes) > 1):
        """in which case randomly select one of these"""
        max_index = np.random.choice(max_indexes)
    else:
        """else just take the maximum value"""
        max_index = max_indexes[0]
    """and return the (random) mode"""
    return max_index

def get_mask_segments_start_and_length(mask):
    k = 0
    splits = []
    is_true = False
    if mask[-1]:
        mask = np.append(mask, False)
    for j in mask:
        if j and not is_true:
            start = k
            is_true = True
        elif not j and is_true:
            length = k - start
            is_true = False
            splits.append((start, length))
        k += 1
    if k == len(mask) and is_true:
        splits.append((start, length))
    return splits

def correct_calipso_phase(phase_in):
    phase_out = np.array(phase_in)
    phase_out = phase_out.astype(float)
    phase_out[phase_out == 1] = 999.
    phase_out[phase_out == 2] = 1.
    phase_out[phase_out == 999.] = 2.
    phase_out[phase_out == 3] = 2.
    phase_out[phase_out == 0.] = 999.
    return phase_out

def correct_cci_phase(phase):
    phase[phase == 0.] = np.nan
    phase[phase > 2.] = np.nan
    return phase

def plot_topography(data, lat, ax):
    """Surface elevation, snow/ice, surface type"""
    topography = data.get('calipsoTop') * 1000.
    topoMax = topography.max()
    ice = data.get('calipsoIce')
    # calipso ice/snow values: 1-100, 101, 102 (not used), 103
    iceNise = data.get('cciNise')
    surface_type = data.get('calipsoTyp')
    ice_masked = ma.masked_inside(ice, 1, 103)
    sea_masked = ma.masked_equal(surface_type, 17)
    topography_masked = topography[ice_masked.mask].copy()
    lat_masked = lat[ice_masked.mask].copy()
    ax2 = ax.twinx()
    ax2.set_ylabel("surface elevation [m]")
    minX = lat.min()
    maxX = lat.max()
    width = (maxX - minX) / len(lat) # bar width
    ax2.set_xlim([minX - width * 0.5, maxX + width * 0.5])
    ax2.set_ylim(0, topoMax * 5)
    if topoMax < 500.:
        yInterval = 100
    else:
        yInterval = 200.
    yticks = np.arange(0, topoMax, yInterval)
    lw = 5.
    ax2.fill_between(lat, topography+10, facecolor="green", alpha=0.5, zorder=1, linewidth=lw, edgecolor="")
    # sea
    if sea_masked.count() < len(sea_masked):
        splits = get_mask_segments_start_and_length(sea_masked.mask)
        for start, length in splits:
            x = lat[start:start+length]
            y = topography[start:start+length]
            ax2.plot(x, y, c="blue", linewidth=lw, zorder=1)
    # ice
    if ice_masked.count() < len(ice_masked):
        splits = get_mask_segments_start_and_length(ice_masked.mask)
        for start, length in splits:
            x = lat[start:start+length]
            y = topography[start:start+length]-10
            ax2.plot(x, y, c="grey", linewidth=3, zorder=1)
    ax2.set_yticks(yticks)

def update_latex_variables(path):
    """loop over all tex files in main folder"""
    for tex in glob.glob(path + "*.tex"):
        """open a test version of each for writing"""
        with open(os.path.splitext(tex)[0] + "_test.tex", 'w') as new_f:
            """open tex file for reading"""
            with open(tex, 'r') as f:
                """loop over all lines of each file"""
                for line in f:
                    """and check whether that line matches the search pattern"""
                    if re.search('insertVariable{.*}[0-9]*\.{1}[0-9]*', line):
                        """if so, loop over all words in line"""
                        for word in line.split():
                            """and check whether that word matches the search pattern"""
                            m = re.search('insertVariable{.*}[0-9]*\.{1}[0-9]*', word)
                            if m:
                                """if so, get the matching word"""
                                found = m.group()
                                """now loop over all latex variables in dictionary"""
                                for variable, value in globals.latex_variables.iteritems():
                                    """and if a key matches the word"""
                                    if variable in found:
                                        value_new = str(value)
                                        """get the old value of the latex variable from the tex file"""
                                        value_old = re.search('[0-9]*\.[0-9]*', found).group()
                                        """replace it with the new value"""
                                        replace = found.replace(value_old, value_new)
                                        """and replace the old word with that new word within the entire line"""
                                        line = line.replace(found, replace)
                    """write each line to the test output file, regardless of whether a match has been found"""
                    new_f.write(line)
        """after having looped over all lines of a file, replace it by its test version"""
        os.remove(tex)
        os.rename(os.path.splitext(tex)[0] + "_test.tex", tex)
