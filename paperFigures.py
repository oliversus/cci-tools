#!/data/osus/Enthought/User/bin/python2.7

from __future__ import division

from analyseCCI import CCI
import numpy as np
from CCITools import buildRGB, plotRGB,\
    plotRGBMulti, greatCircle, collocateCciAndCalipso, \
    plotCciCalipsoCollocation, plotCCI, plotCCIMulti, \
    update_latex_variables, calculate_statistics
import sys
import numpy.ma as ma
from pyhdf.SD import SD, SDC
from sys import argv
import math
import globals
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

"""initialise global variables"""
globals.init()

if len(argv) > 1:
    delLon = argv[1]
    delLat = argv[2]
    if argv[3] == "True":
        doRGB = True
    elif argv[3] == "False":
        doRGB = False
    else:
        print "ERROR: 3rd argument doRGB should be [True/False]."
        sys.exit()
    sceneTime = argv[4] # 1) 07221915 2) 07270810 3) 07230021 4) 07222058
    if sceneTime not in globals.sceneTimes:
    #if argv[4] != '07221915' and argv[4] != '07270810' and argv[4] != '07230021' and argv[4] != '07222058':
        print "ERROR: choose correct study date ('07221915' or '07270810' or '07230021' or '07222058')"
        sys.exit()
    if argv[5] == "True":
        corrected = True
    elif argv[5] == "False":
        corrected = False
    else:
        print "ERROR: 5th argument corrected should be [True/False]."
        sys.exit()
    if argv[6] == "True":
        plotCot = True
    elif argv[6] == "False":
        plotCot = False
    else:
        print "ERROR: 6th argument plotCot should be [True/False]."
        sys.exit()
    if argv[7] == "True":
        plotCalipso = True
    elif argv[7] == "False":
        plotCalipso = False
    else:
        print "ERROR: 7th argument plotCot should be [True/False]."
        sys.exit()
else:
    delLat = "0.1"
    delLon = "0.1"
    doRGB = False
    sceneTime = '07222058'
    corrected = False
    plotCot = False
    plotCalipso = True

plot_variables = False
plot_statistics = False
globals.sceneTime = sceneTime
delLat_delLon = delLat + "_" + delLon
month = sceneTime[0:2]
day = sceneTime[2:4]
hour = sceneTime[4:6]
minute = sceneTime[6:8]

calipsoPath1km = globals.data_folder + "calipso_1km_" + sceneTime + ".hdf"
calipsoPath5km = globals.data_folder + "calipso_5km_" + sceneTime + ".hdf"

hdf = SD(calipsoPath5km, SDC.READ)
# Read geolocation
lat = hdf.select('Latitude')
lon = hdf.select('Longitude')
cod = hdf.select('Column_Optical_Depth_Cloud_532')
codCum = hdf.select('Feature_Optical_Depth_532')
ctp = hdf.select('Layer_Top_Pressure')
ctpBot = hdf.select('Layer_Base_Pressure')
ctt = hdf.select('Layer_Top_Temperature')
cth = hdf.select('Layer_Top_Altitude')
cloudFlag = hdf.select('Feature_Classification_Flags')
sfcIce = hdf.select('NSIDC_Surface_Type')
sfcElev = hdf.select('DEM_Surface_Elevation')
sfcType = hdf.select('IGBP_Surface_Type')
calipsoLat = lat[:,1] # 4208,
calipsoLon = lon[:,1] # 4208,
calipsoCOD = cod[:,0] # 4208,
calipsoCTP = ctp[:,:] # 4208, 10
calipsoCTPBot = ctpBot[:,:] # 4208, 10
calipsoCTT = ctt[:,:] # 4208, 10
calipsoCTH = cth[:,:] # 4208, 10
calipsoFCF = cloudFlag[:,:] # 4208, 10: 0 is top, 9 bottom layer
calipsoICE = sfcIce[:,:]
calipsoTOP = sfcElev[:,:]
calipsoTYP = sfcType[:,:]
calipsoCODLayer = codCum[:,:] # 4208, 10: 0 is top, 9 bottom layer
calipsoData = {'lat': calipsoLat, 'lon': calipsoLon,
               'cod': calipsoCOD, 'codLayered': calipsoCODLayer,
               'ctp': calipsoCTP, 'ctpBot': calipsoCTPBot,
               'ctt': calipsoCTT, 'cth': calipsoCTH,
               'fcf': calipsoFCF, 'ice': calipsoICE,
               'top': calipsoTOP, 'typ': calipsoTYP}

delLatStr = str(delLat); delLatStr = delLatStr.replace(".", "")
delLonStr = str(delLon); delLonStr = delLonStr.replace(".", "")

l2_primary_prefix = "cci_l2_primary_"
N18PrimaryResampledName = l2_primary_prefix + "n18_" + sceneTime + "_resampled_" + str(delLon) + "_" + str(delLat) + ".nc"
N18SecondaryResampledName = N18PrimaryResampledName.replace("primary", "secondary")
MYDPrimaryResampledName = N18PrimaryResampledName.replace("n18", "myd")
MYDSecondaryResampledName = MYDPrimaryResampledName.replace("primary", "secondary")
ENVPrimaryResampledName = N18PrimaryResampledName.replace("n18", "env")
ENVSecondaryResampledName = ENVPrimaryResampledName.replace("primary", "secondary")

# NOAA18 paths and data
print "Reading NOAA18 data"
pathL2PriN18 = globals.data_folder + l2_primary_prefix + "n18_" + sceneTime + ".nc"
pathL2SecN18 = pathL2PriN18.replace("primary", "secondary")
#priN18 = CCI(pathL2PriN18)
#secN18 = CCI(pathL2SecN18)

# MODIS AQUA paths and data
print "Reading MODIS AQUA data"
pathL2PriMYD = globals.data_folder + l2_primary_prefix + "myd_" + sceneTime + ".nc"
pathL2SecMYD = pathL2PriMYD.replace("primary", "secondary")
#priMYD = CCI(pathL2PriMYD)
#secMYD = CCI(pathL2SecMYD)

# ENVISAT AATSR paths and data
print "Reading ENVISAT AATSR data"
pathL2PriENV = globals.data_folder + l2_primary_prefix + "env_" + sceneTime + ".nc"
pathL2SecENV = pathL2PriENV.replace("primary", "secondary")
#priENV = CCI(pathL2PriENV)
#secENV = CCI(pathL2SecENV)

# N18 paths and data, RESAMPLED
print "Reading resampled N18 data"
pathL2N18PrimaryResampled = globals.data_folder + N18PrimaryResampledName
N18PrimaryResampled = CCI(pathL2N18PrimaryResampled)
N18Masked = CCI(pathL2N18PrimaryResampled)
pathL2N18SecondaryResampled = globals.data_folder + N18SecondaryResampledName
N18SecondaryResampled = CCI(pathL2N18SecondaryResampled)

# MODIS AQUA paths and data, RESAMPLED
print "Reading resampled MODIS AQUA data"
pathL2MYDPrimaryResampled = globals.data_folder + MYDPrimaryResampledName
MYDPrimaryResampled = CCI(pathL2MYDPrimaryResampled)
pathL2MYDSecondaryResampled = globals.data_folder + MYDSecondaryResampledName
MYDSecondaryResampled = CCI(pathL2MYDSecondaryResampled)

# ENVISAT AATSR paths and data, RESAMPLED
print "Reading resampled ENVISAT AATSR data"
pathL2ENVPrimaryResampled = globals.data_folder + ENVPrimaryResampledName
ENVPrimaryResampled = CCI(pathL2ENVPrimaryResampled)
pathL2ENVSecondaryResampled = globals.data_folder + ENVSecondaryResampledName
ENVSecondaryResampled = CCI(pathL2ENVSecondaryResampled)

# subset borders in lat/lon
if sceneTime == '07221915':
    centrePoint = [72.5, -111.]
    boundingBox = [-179., 0., 40, 90]
    poly_lats = [58.5, 78, 83.5, 63]
    poly_lons = [-121, -73, -63., -128.5]
elif sceneTime == '07270810':
    centrePoint = [73., 55.]
    boundingBox = [-10., 130., 45., 90.]
    poly_lats = [62.5, 78., 83.5, 63.5]
    poly_lons = [43, 105, 110, 37.5]
elif sceneTime == '07230021':
    centrePoint = [71., 173.]
    boundingBox = [140., -160., 45., 90.]
    poly_lats = [63, 79., 85, 65]
    poly_lons = [45, 90, 70, 30]
elif sceneTime == '07222058':
    centrePoint = [62., -125.]
    boundingBox = [-180., -100., 45., 90.]
    poly_lats = [63, 72.5, 85, 65]
    poly_lons = [45, 90, 70, 30]

# get all variables
#MYDSlice = priMYD.getAllVariables(doSlice=True, boundingBox=boundingBox)
#secMYD.getAllVariables(doSlice=True, boundingBox=boundingBox, primary=False, boxSlice=MYDSlice)
print "Getting all variables: N18 resampled"
N18PrimaryResampled.getAllVariables()
N18Masked.getAllVariables()
N18SecondaryResampled.getAllVariables()
print "Getting all variables: MODIS resampled"
MYDPrimaryResampled.getAllVariables()
MYDSecondaryResampled.getAllVariables()
print "Getting all variables: Envisat resampled"
ENVPrimaryResampled.getAllVariables()
ENVSecondaryResampled.getAllVariables()

# mask all resampled pixels with cc_total < 1 to exclude fractional cloud coverage
N18ResampledCloudMask = ma.masked_less(N18PrimaryResampled.cc_total, 1.).mask
MYDResampledCloudMask = ma.masked_less(MYDPrimaryResampled.cc_total, 1.).mask
ENVResampledCloudMask = ma.masked_less(ENVPrimaryResampled.cc_total, 1.).mask
CloudMask = N18ResampledCloudMask + MYDResampledCloudMask + ENVResampledCloudMask

# build mask of all pixels out of study area, i.e. where any sensor has no reflectance data
N18ReflMask = N18SecondaryResampled.reflectance_in_channel_no_1.mask
MYDReflMask = MYDSecondaryResampled.reflectance_in_channel_no_1.mask
ENVReflMask = ENVSecondaryResampled.reflectance_in_channel_no_1.mask
ReflMask = N18ReflMask + MYDReflMask + ENVReflMask

SuperMask = CloudMask + ReflMask
#ReflMask = SuperMask

N18Masked.maskAllVariables(ReflMask)
poly_lats = [62.5, 78., 83.5, 63.5]
poly_lats[0] = N18Masked.lat.min()
#poly_lats[1] = min(N18Masked.lat[N18Masked.lon == N18Masked.lon.max()])
poly_lats[2] = N18Masked.lat.max()
poly_lats[3] = max(N18Masked.lat[N18Masked.lon == N18Masked.lon.min()])
poly_lons = [43, 105, 110, 37.5]
poly_lons[0] = min(N18Masked.lon[N18Masked.lat == N18Masked.lat.min()])
poly_lons[1] = N18Masked.lon.max()
poly_lons[2] = max(N18Masked.lon[N18Masked.lat == N18Masked.lat.max()])
poly_lons[3] = N18Masked.lon.min()

if sceneTime == '07222058':
    poly_lats[1] = 72.5
    poly_lons[1] = -101
elif sceneTime == '07221915':
    poly_lats[1] = 78
    poly_lons[1] = -73
elif sceneTime == '07270810':
    poly_lats[0] = 60
    poly_lons[0] = 45
    poly_lats[1] = 78
    poly_lons[1] = 95.5

#################################################################################################
# plot variables and calculate statistics
# 1) compare CCI products from N18, MYD, and AATSR

    # a) RGB, ideally highly resolved MODIS with 0.6, 0.8, and 1.6 (1.6 has missing scan lines though)
if doRGB:    
    print "RGB: started."
    platform = "MYD"
    colourTupleMYD = buildRGB(MYDPrimaryResampled, MYDSecondaryResampled, platform)
    RGBName = globals.figuresDir + "RGB_" + platform + "_" + delLatStr + "x" + delLonStr + "_" + sceneTime + ".png"
    plotRGB(RGBName, colourTupleMYD, MYDSecondaryResampled.lat, MYDSecondaryResampled.lon, MYDSecondaryResampled.reflectance_in_channel_no_1,
            centrePoint[0], centrePoint[1])
    platform = "N18"
    colourTupleN18 = buildRGB(N18PrimaryResampled, N18SecondaryResampled, platform)
    RGBName = globals.figuresDir + "RGB_" + platform + "_" + delLatStr + "x" + delLonStr + "_" + sceneTime + ".png"
    plotRGB(RGBName, colourTupleN18, N18SecondaryResampled.lat, N18SecondaryResampled.lon, N18SecondaryResampled.reflectance_in_channel_no_1,
            centrePoint[0], centrePoint[1])
    platform = "ENV"
    colourTupleENV = buildRGB(ENVPrimaryResampled, ENVSecondaryResampled, platform)
    RGBName = globals.figuresDir + "RGB_" + platform + "_" + delLatStr + "x" + delLonStr + "_" + sceneTime + ".png"
    plotRGB(RGBName, colourTupleENV, ENVSecondaryResampled.lat, ENVSecondaryResampled.lon, ENVSecondaryResampled.reflectance_in_channel_no_1,
            centrePoint[0], centrePoint[1])
    colourTupleMulti = np.concatenate((colourTupleN18[..., np.newaxis], colourTupleMYD[..., np.newaxis], colourTupleENV[..., np.newaxis]), axis=2)
    RGBName = globals.figuresDir + "RGB_multi_" + delLatStr + "x" + delLonStr + "_" + sceneTime + ".png"
    plotRGBMulti(RGBName, colourTupleMulti, N18SecondaryResampled.lat, N18SecondaryResampled.lon, poly_lats, poly_lons, N18SecondaryResampled.reflectance_in_channel_no_1,
                calipsoLat, calipsoLon,
                centrePoint[0], centrePoint[1])    
    print "RGB: done."    

print "Plotting data."

# 1a) RGB, ideally highly resolved MODIS with 0.6, 0.8, and 1.6 (1.6 has missing scan lines though)

if plot_variables:
    # 1b) calculate min/mean/max time difference between grid boxes = TIMEDIFF (no plot)
    variable = "ctp"
    plotCCI(N18PrimaryResampled, MYDPrimaryResampled, ENVPrimaryResampled, boundingBox, centrePoint, variable,
            'N18', colourMin=0, colourMax=1000)
    plotCCI(N18PrimaryResampled, MYDPrimaryResampled, ENVPrimaryResampled, boundingBox, centrePoint, variable,
            'MYD', colourMin=0, colourMax=1000)
    plotCCI(N18PrimaryResampled, MYDPrimaryResampled, ENVPrimaryResampled, boundingBox, centrePoint, variable,
            'ENV', colourMin=0, colourMax=1000)
    plotCCIMulti(N18PrimaryResampled, MYDPrimaryResampled, ENVPrimaryResampled, boundingBox, centrePoint, variable,
                 poly_lats, poly_lons, colourMin=0, colourMax=1000)

    variable = "cot"
    plotCCI(N18PrimaryResampled, MYDPrimaryResampled, ENVPrimaryResampled, boundingBox, centrePoint, variable,
            'N18', colourMin=0, colourMax=50)
    plotCCI(N18PrimaryResampled, MYDPrimaryResampled, ENVPrimaryResampled, boundingBox, centrePoint, variable,
            'MYD', colourMin=0, colourMax=50)
    plotCCI(N18PrimaryResampled, MYDPrimaryResampled, ENVPrimaryResampled, boundingBox, centrePoint, variable,
            'ENV', colourMin=0, colourMax=50)
    plotCCIMulti(N18PrimaryResampled, MYDPrimaryResampled, ENVPrimaryResampled, boundingBox, centrePoint, variable,
                 poly_lats, poly_lons, colourMin=0, colourMax=50)

    variable = "cer"
    plotCCI(N18PrimaryResampled, MYDPrimaryResampled, ENVPrimaryResampled, boundingBox, centrePoint, variable,
            'N18', colourMin=0, colourMax=50)
    plotCCI(N18PrimaryResampled, MYDPrimaryResampled, ENVPrimaryResampled, boundingBox, centrePoint, variable,
            'MYD', colourMin=0, colourMax=50)
    plotCCI(N18PrimaryResampled, MYDPrimaryResampled, ENVPrimaryResampled, boundingBox, centrePoint, variable,
            'ENV', colourMin=0, colourMax=50)
    plotCCIMulti(N18PrimaryResampled, MYDPrimaryResampled, ENVPrimaryResampled, boundingBox, centrePoint, variable,
                 poly_lats, poly_lons, colourMin=0, colourMax=50)

    # plot uncertainties
    figure_name = globals.figuresDir + sceneTime + "_uncertainties_absolute.png" # "_uncertainties_percent.png"
    fig1 = plt.figure(figsize=(15, 15))
    gs1 = gridspec.GridSpec(2, 2)
    gs1.update(wspace=0.15, hspace=0.15) # set the spacing between axes.
    variable = 'ctp'
    variable_unc = variable + '_uncertainty'
    ax = plt.subplot(gs1[0])
    ax.title.set_text('CTP')
    input = getattr(MYDPrimaryResampled, variable_unc) # 100. * getattr(MYDPrimaryResampled, variable_unc) / getattr(MYDPrimaryResampled, variable)
    if sceneTime == globals.NA2:
        globals.latex_variables["NA2_ctp_unc_lt10"] = 100. * np.round(np.sum(input < 10.) / ma.count(input), 1)
        foo = input[~input.mask]
        globals.latex_variables["NA2_ctp_unc_median"] = np.round(np.nanmedian(foo), 1)
    plotCCI(N18PrimaryResampled, MYDPrimaryResampled, ENVPrimaryResampled, boundingBox, centrePoint, variable,
            'MYD', input=input, colourMin=0, colourMax=50, create_figure=False)
    variable = 'cot'
    variable_unc = variable + '_uncertainty'
    ax = plt.subplot(gs1[1])
    #ax = fig1.add_subplot(2, 2, 2)
    ax.title.set_text('COT')
    input = getattr(MYDPrimaryResampled, variable_unc) # 100. * getattr(MYDPrimaryResampled, variable_unc) / getattr(MYDPrimaryResampled, variable)
    if sceneTime == globals.NA2:
        foo = input[~input.mask]
        globals.latex_variables["NA2_cot_unc_median"] = np.round(np.nanmean(foo), 1)
    plotCCI(N18PrimaryResampled, MYDPrimaryResampled, ENVPrimaryResampled, boundingBox, centrePoint, variable,
            'MYD', input=input, colourMin=0, colourMax=50, create_figure=False)
    variable = 'cer'
    variable_unc = variable + '_uncertainty'
    ax = plt.subplot(gs1[2])
    ax.title.set_text('CER')
    input = getattr(MYDPrimaryResampled, variable_unc) # 100. * getattr(MYDPrimaryResampled, variable_unc) / getattr(MYDPrimaryResampled, variable)
    if sceneTime == globals.NA2:
        foo = input[~input.mask]
        globals.latex_variables["NA2_cer_unc_median"] = np.round(np.nanmedian(foo), 1)
    plotCCI(N18PrimaryResampled, MYDPrimaryResampled, ENVPrimaryResampled, boundingBox, centrePoint, variable,
            'MYD', input=input, colourMin=0, colourMax=50, create_figure=False)
    variable = 'cc_total_unc'
    ax = plt.subplot(gs1[3])
    ax.title.set_text('Cloud mask')
    input = getattr(MYDPrimaryResampled, variable)
    if sceneTime == globals.NA2:
        foo = input[~input.mask]
        globals.latex_variables["NA2_cmask_unc_median"] = np.round(np.nanmedian(foo), 1)
    plotCCI(N18PrimaryResampled, MYDPrimaryResampled, ENVPrimaryResampled, boundingBox, centrePoint, variable,
            'MYD', colourMin=0, colourMax=50, create_figure=False)
    plt.savefig(figure_name, bbox_inches='tight')

########################################
""" Calipso collocation with CCI grid"""

"""mask all CCI pixels that do not have reflectance values for all 3 sensors"""
N18PrimaryResampled.maskAllVariables(ReflMask)
MYDPrimaryResampled.maskAllVariables(ReflMask)
ENVPrimaryResampled.maskAllVariables(ReflMask)
"""the maximum distance between Calipso and a grid box center is given by
    the great circle distance between the grid box center and one of its corners
    which is half of the grid box diagonal, which can be calculated with Pythagoras"""
"""The lat/lon spacing"""
dello = float(delLon)
della = float(delLat)
"""and the grid box longitude/latitude width"""
boxLonWidth = dello * greatCircle(centrePoint[1], centrePoint[0], centrePoint[1] + 1., centrePoint[0])
boxLatWidth = della * greatCircle(centrePoint[1], centrePoint[0], centrePoint[1], centrePoint[0] + 1.)
"""give the grid box diagonal"""
boxDiag = (boxLonWidth**2 + boxLatWidth**2)**0.5
"""half of which is the maximum distance of Calipso to any grid box center"""
maxDistance = boxDiag / 2.
"""The maximum distance is input to the collocation method,
    returning collocated CCI and Calipso data."""
collocateN18 = collocateCciAndCalipso(N18PrimaryResampled, calipsoData, maxDistance, corrected)
collocateMYD = collocateCciAndCalipso(MYDPrimaryResampled, calipsoData, maxDistance, corrected)
collocateENV = collocateCciAndCalipso(ENVPrimaryResampled, calipsoData, maxDistance, corrected)

"""Plot collocated data for COT, CTP, and CTT."""
figurePath = globals.figuresDir + "calipsoVsCci_" + sceneTime
if plotCot:
    figurePath += "_cot"
else:
    figurePath += "_nocot"
if corrected:
    figurePath += "_correctedCtp"
else:
    figurePath += "_uncorrectedCtp"
figurePath += ".png"
if plotCalipso:
    plotCciCalipsoCollocation(collocateN18, collocateMYD, collocateENV, figurePath, sceneTime, plotCot)

    # b) calculate min/mean/max time difference between grid boxes = TIMEDIFF (no plot)
timeDiff = ((MYDPrimaryResampled.time - math.floor(MYDPrimaryResampled.time.min())) \
            - (N18PrimaryResampled.time - math.floor(N18PrimaryResampled.time.min()))) * 24 * 60
if True:
    print "The minimum/mean/maximum time difference between MODIS AQUA and NOAA18 is: ", \
        round(timeDiff.min(), 2), "min/", round(np.mean(timeDiff), 2), "min/", round(timeDiff.max(), 2), "min"

    # c) calculate min/mean/max solzen and satzen for all 3 sensors
if True:
    print "Min solzen MYD  = ", round(MYDPrimaryResampled.solar_zenith_view_no1.min(), 2)
    print "Mean solzen MYD = ", round(np.mean(MYDPrimaryResampled.solar_zenith_view_no1), 2)
    print "Max solzen MYD  = ", round(MYDPrimaryResampled.solar_zenith_view_no1.max(), 2)
        # check that all sensors have same amount of data, otherwise create mask --> print n values

    # d) plot ctp and cot for N18, MYD, AATSR

    # e) plot histogram with normal (?) distribution fit

    # f) calculate statistical moments: mean, median, std, skew, kurtosis
         # plot boxplots?

    # g) calculate significance of differences in variance and mean

    # h) calculate residuals and analyse dependence on third variable
        # RES = CTP_MYD - CTP_N18 -> mean/median/std
        # linear regression of RES with COT or TIMEDIFF
            # scatter plot RES vs COT/TIMEDIFF?

#sys.exit()


# calculate statistics for study area polygon only

if plot_statistics:

    ReflMask = SuperMask
    N18PrimaryResampled.maskAllVariables(ReflMask)
    N18SecondaryResampled.maskAllVariables(ReflMask)
    MYDPrimaryResampled.maskAllVariables(ReflMask)
    MYDSecondaryResampled.maskAllVariables(ReflMask)
    ENVPrimaryResampled.maskAllVariables(ReflMask)
    ENVSecondaryResampled.maskAllVariables(ReflMask)

    nbins = 30
    ttest_threshold = 0.01
    variable = 'ctp'
    x = getattr(N18PrimaryResampled, variable).ravel().compressed()
    y = getattr(MYDPrimaryResampled, variable).ravel().compressed()
    z = getattr(ENVPrimaryResampled, variable).ravel().compressed()
    if stats.ttest_ind(x, y).pvalue > ttest_threshold:
        print variable + " ttest for N18 + MYD is > " + str(ttest_threshold) + ": " + str(stats.ttest_ind(x, y).pvalue)
    if stats.ttest_ind(x, z).pvalue > ttest_threshold:
        print variable + " ttest for N18 + ENV is > " + str(ttest_threshold) + ": " + str(stats.ttest_ind(x, z).pvalue)
    if stats.ttest_ind(y, z).pvalue > ttest_threshold:
        print variable + " ttest for MYD + ENV is > " + str(ttest_threshold) + ": " + str(stats.ttest_ind(y, z).pvalue)
    if sceneTime == globals.NA2:
        calculate_statistics(x, variable, 'N18')
        calculate_statistics(y, variable, 'MYD')
        calculate_statistics(z, variable, 'ENV')
    data = np.vstack([x,y,z]).T
    bins = np.linspace(200, 1100, nbins)
    fig1 = plt.figure(figsize=(30, 10))
    ax=fig1.add_subplot(1,3,1)
    plt.hist(data, range=(0, 1100), normed=True, label=['N18','MYD','ENV'], bins=bins, color=['r', 'aqua', 'orange'])
    plt.legend(loc=2, title="(a)")
    ax.set_xlabel(variable.upper())
    variable = 'cot'
    x = getattr(N18PrimaryResampled, variable).ravel().compressed()
    y = getattr(MYDPrimaryResampled, variable).ravel().compressed()
    z = getattr(ENVPrimaryResampled, variable).ravel().compressed()
    if stats.ttest_ind(x, y).pvalue > ttest_threshold:
        print variable + " ttest for N18 + MYD is > " + str(ttest_threshold) + ": " + str(stats.ttest_ind(x, y).pvalue)
    if stats.ttest_ind(x, z).pvalue > ttest_threshold:
        print variable + " ttest for N18 + ENV is > " + str(ttest_threshold) + ": " + str(stats.ttest_ind(x, z).pvalue)
    if stats.ttest_ind(y, z).pvalue > ttest_threshold:
        print variable + " ttest for MYD + ENV is > " + str(ttest_threshold) + ": " + str(stats.ttest_ind(y, z).pvalue)
    if sceneTime == globals.NA2:
        calculate_statistics(x, variable, 'N18')
        calculate_statistics(y, variable, 'MYD')
        calculate_statistics(z, variable, 'ENV')
    data = np.vstack([x,y,z]).T
    bins = np.linspace(0, 50, nbins)
    ax=fig1.add_subplot(1,3,2)
    plt.hist(data, range=(0, 50), normed=True, label=['N18','MYD','ENV'], bins=bins, color=['r', 'aqua', 'orange'])
    plt.legend(loc=1, title="(b)")
    ax.set_xlabel(variable.upper())
    variable = 'cer'
    x = getattr(N18PrimaryResampled, variable).ravel().compressed()
    y = getattr(MYDPrimaryResampled, variable).ravel().compressed()
    z = getattr(ENVPrimaryResampled, variable).ravel().compressed()
    if stats.ttest_ind(x, y).pvalue > ttest_threshold:
        print variable + " ttest for N18 + MYD is > " + str(ttest_threshold) + ": " + str(stats.ttest_ind(x, y).pvalue)
    if stats.ttest_ind(x, z).pvalue > ttest_threshold:
        print variable + " ttest for N18 + ENV is > " + str(ttest_threshold) + ": " + str(stats.ttest_ind(x, z).pvalue)
    if stats.ttest_ind(y, z).pvalue > ttest_threshold:
        print variable + " ttest for MYD + ENV is > " + str(ttest_threshold) + ": " + str(stats.ttest_ind(y, z).pvalue)
    if sceneTime == globals.NA2:
        calculate_statistics(x, variable, 'N18')
        calculate_statistics(y, variable, 'MYD')
        calculate_statistics(z, variable, 'ENV')
    data = np.vstack([x,y,z]).T
    bins = np.linspace(0, 50, nbins)
    ax=fig1.add_subplot(1,3,3)
    plt.hist(data, range=(0, 50), normed=True, label=['N18','MYD','ENV'], bins=bins, color=['r', 'aqua', 'orange'])
    plt.legend(loc=1, title="(c)")
    ax.set_xlabel(variable.upper())

    # differences between retrievals: N18 - MYD, N18 - ENV, MYD - ENV
    variable = 'ctp'
    x = getattr(N18PrimaryResampled, variable).ravel().compressed()
    y = getattr(MYDPrimaryResampled, variable).ravel().compressed()
    z = getattr(ENVPrimaryResampled, variable).ravel().compressed()
    a = x - y
    b = x - z
    c = y - z
    if sceneTime == globals.NA2:
        calculate_statistics(a, variable + 'd', 'N18')
        calculate_statistics(b, variable + 'd', 'MYD')
        calculate_statistics(c, variable + 'd', 'ENV')
    # data = np.vstack([a, b, c]).T
    # bins = np.linspace(-100, 100, nbins)
    # fig1 = plt.figure(figsize = (10, 10))
    # ax=fig1.add_subplot(2,3,4)
    # plt.hist(data, range=(0, 1100), normed=True, label=['N18 - MYD','N18 - ENV','MYD - ENV'], bins=bins, color=['r', 'aqua', 'orange'])
    # plt.legend(loc=2, title="(b)")
    # ax.set_xlabel("$\Delta$ " + variable.upper())
    variable = 'cot'
    x = getattr(N18PrimaryResampled, variable).ravel().compressed()
    y = getattr(MYDPrimaryResampled, variable).ravel().compressed()
    z = getattr(ENVPrimaryResampled, variable).ravel().compressed()
    a = x - y
    b = x - z
    c = y - z
    if sceneTime == globals.NA2:
        calculate_statistics(a, variable + 'd', 'N18')
        calculate_statistics(b, variable + 'd', 'MYD')
        calculate_statistics(c, variable + 'd', 'ENV')
    # data = np.vstack([a, b, c]).T
    # bins = np.linspace(-10, 10, nbins)
    # ax=fig1.add_subplot(2,3,5)
    # plt.hist(data, range=(0, 50), normed=True, label=['N18 - MYD','N18 - ENV','MYD - ENV'], bins=bins, color=['r', 'aqua', 'orange'])
    # ax.set_xlabel("$\Delta$ " + variable.upper())
    # plt.legend(loc=1, title="(d)")
    variable = 'cer'
    x = getattr(N18PrimaryResampled, variable).ravel().compressed()
    y = getattr(MYDPrimaryResampled, variable).ravel().compressed()
    z = getattr(ENVPrimaryResampled, variable).ravel().compressed()
    a = x - y
    b = x - z
    c = y - z
    if sceneTime == globals.NA2:
        calculate_statistics(a, variable + 'd', 'N18')
        calculate_statistics(b, variable + 'd', 'MYD')
        calculate_statistics(c, variable + 'd', 'ENV')
    # data = np.vstack([a, b, c]).T
    # bins = np.linspace(-10, 10, nbins)
    # ax=fig1.add_subplot(2,3,6)
    # plt.hist(data, range=(0, 50), normed=True, label=['N18 - MYD','N18 - ENV','MYD - ENV'], bins=bins, color=['r', 'aqua', 'orange'])
    # ax.set_xlabel("$\Delta$ " + variable.upper())
    # plt.legend(loc=1, title="(f)")
    plt.savefig(globals.figuresDir + sceneTime + '_histograms.png', bbox_inches='tight')


    # ax=fig1.add_subplot(2,1,2)
    # #plt.hist(data[:,0], range=(0, 1100), normed=True, bins=bins)
    # plt.hist(x, range=(0, 1100), normed=True, bins=bins)

    input = abs(getattr(MYDPrimaryResampled, variable) - getattr(N18PrimaryResampled, variable))
    #input.mask = input.mask + ma.masked_greater(input, input.mean() + input.std() * 2).mask

    # plotCCI(N18PrimaryResampled, MYDPrimaryResampled, boundingBox, centrePoint, variable,
    #         'MYD', mask = AllSensorsMaskCombined)
    # print "MYD " + variable + " mean = " + str(round(getattr(MYDPrimaryResampled, variable).mean(), 2))
    # print "MYD " + variable + " stdev = " + str(round(getattr(MYDPrimaryResampled, variable).std(), 2))
    # input = abs(getattr(MYDPrimaryResampled, variable) - getattr(N18PrimaryResampled, variable))
    # print [round(x, 2) for x in minMax(input)]
    # plotCCI(N18PrimaryResampled, MYDPrimaryResampled, boundingBox, centrePoint, variable,
    #         'MYDResampledMinusNOAA18', input = input, mask = input.mask)
    # input.mask = input.mask + ma.masked_greater(input, input.mean() + input.std() * 2).mask
    # plotCCI(N18PrimaryResampled, MYDPrimaryResampled, boundingBox, centrePoint, variable,
    #         'MYDResampledMinusNOAA18', input = input, mask = input.mask)
    # print "MYD-N18 " + variable + " mean = " + str(round(input.mean(), 2))
    # print "MYD-N18 " + variable + " stdev = " + str(round(input.std(), 2))
    # plt.scatter(MYDPrimaryResampled.ctp, input)
    # plt.draw()
    # plt.scatter(timeDiff, input)
    # plt.draw()

    # ttests
    # N18 and MYD
    print stats.ttest_ind(getattr(N18PrimaryResampled, variable).ravel(), getattr(MYDPrimaryResampled, variable).ravel())
    # N18 and ENV
    print stats.ttest_ind(getattr(N18PrimaryResampled, variable).ravel(), getattr(ENVPrimaryResampled, variable).ravel())
    # MYD and ENV
    print stats.ttest_ind(getattr(MYDPrimaryResampled, variable).ravel(), getattr(ENVPrimaryResampled, variable).ravel())

print "updating latex variables"
print globals.latex_variables
path = globals.main_folder
update_latex_variables(path)
print "...done"
