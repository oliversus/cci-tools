#!/data/osus/Enthought/User/bin/python2.7

from analyseCCI import CCI
import numpy as np
from CCITools import buildRGB, plotRGB,\
    plotRGBMulti, greatCircle, collocateCciAndCalipso, \
    plotCciCalipsoCollocation, plotCCI, minMax
import sys
import numpy.ma as ma
from pyhdf.SD import SD, SDC
from sys import argv
import math
import matplotlib.pyplot as plt

if len(argv) > 1:
    delLon = argv[1]
    delLat = argv[2]
    doRGB  = argv[3]
    if argv[3] == "True":
        doRGB = True
    elif argv[3] == "False":
        doRGB = False
    else:
        print "ERROR: 3rd argument should be [True/False]."
        sys.exit()
else:
    delLat = "0.1"
    delLon = "0.1"
    doRGB = False


calipsoPath1km = "/cmsaf/cmsaf-cld1/thanschm/VALIDATION/DATASETS/CALIOP/1kmClay/2008/07/22/CAL_LID_L2_01kmCLay-ValStage1-V3-01.2008-07-22T18-41-41ZD.hdf"
calipsoPath5km = "/cmsaf/cmsaf-cld1/thanschm/VALIDATION/DATASETS/CALIOP/2008/07/22/CAL_LID_L2_05kmCLay-Prov-V3-01.2008-07-22T18-41-41ZD.hdf"
hdf = SD(calipsoPath5km, SDC.READ)
# Read geolocation
lat = hdf.select('Latitude')
lon = hdf.select('Longitude')
cod = hdf.select('Column_Optical_Depth_Cloud_532')
codCum = hdf.select('Feature_Optical_Depth_532')
ctp = hdf.select('Layer_Top_Pressure')
ctt = hdf.select('Layer_Top_Temperature')
cloudFlag = hdf.select('Feature_Classification_Flags')
sfcIce = hdf.select('NSIDC_Surface_Type')
sfcElev = hdf.select('DEM_Surface_Elevation')
sfcType = hdf.select('IGBP_Surface_Type')
calipsoLat = lat[:,1] # 4208,
calipsoLon = lon[:,1] # 4208,
calipsoCOD = cod[:,0] # 4208,
calipsoCTP = ctp[:,:] # 4208, 10
calipsoCTT = ctt[:,:] # 4208, 10
calipsoFCF = cloudFlag[:,:] # 4208, 10: 0 is top, 9 bottom layer
calipsoICE = sfcIce[:,:]
calipsoTOP = sfcElev[:,:]
calipsoTYP = sfcType[:,:]
calipsoCODLayer = codCum[:,:] # 4208, 10: 0 is top, 9 bottom layer
calipsoData = {'lat': calipsoLat, 'lon': calipsoLon,
               'cod': calipsoCOD, 'codLayered': calipsoCODLayer,
               'ctp': calipsoCTP, 'ctt': calipsoCTT,
               'fcf': calipsoFCF, 'ice': calipsoICE,
               'top': calipsoTOP, 'typ': calipsoTYP}

figuresDir = "/cmsaf/esa_doku/ESA_Cloud_cci/publications/CC4CL_paper/figures/"
mainL1 = "/cmsaf/cmsaf-cld7/esa_cloud_cci/data/v2.0/L1/"
mainL2 = "/cmsaf/cmsaf-cld7/esa_cloud_cci/data/v2.0/L2/"
delLatStr = str(delLat); delLatStr = delLatStr.replace(".", "")
delLonStr = str(delLon); delLonStr = delLonStr.replace(".", "")
N18PrimaryResampledName = "N18_Resampled_2008-07-22-1851_lat" + delLat + "lon" + delLon + "_primary.nc"
MYDPrimaryResampledName = "MYD_Resampled_2008-07-22-1915_lat" + delLat + "lon" + delLon + "_primary.nc"
ENVPrimaryResampledName = "ENV_Resampled_2008-07-22-1844_lat" + delLat + "lon" + delLon + "_primary.nc"
N18SecondaryResampledName = "N18_Resampled_2008-07-22-1851_lat" + delLat + "lon" + delLon + "_secondary.nc"
MYDSecondaryResampledName = "MYD_Resampled_2008-07-22-1915_lat" + delLat + "lon" + delLon + "_secondary.nc"
ENVSecondaryResampledName = "ENV_Resampled_2008-07-22-1844_lat" + delLat + "lon" + delLon + "_secondary.nc"

# NOAA18 paths and data
print "Reading NOAA18 data"
pathL2PriN18 = mainL2 + "20080722185100-ESACCI-L2_CLOUD-CLD_PRODUCTS-AVHRRGAC-NOAA18-fv2.0.nc"
pathL2SecN18 = mainL2 + "ECC_GAC_avhrr_noaa18_99999_20080722T1851289Z_20080722T2046134Z.secondary.nc"
priN18 = CCI(pathL2PriN18)
secN18 = CCI(pathL2SecN18)

# MODIS AQUA paths and data
print "Reading MODIS AQUA data"
pathL2PriMYD = mainL2 + "MYD_merged_20080722_19151920_primary.nc" #"MYD20080722_1915.nc"
pathL2SecMYD = mainL2 + "MYD_merged_20080722_19151920_secondary.nc" #"MYD021KM.A2008204.1915.006.2012069115248.bspscs_000500694537.secondary.nc"
priMYD = CCI(pathL2PriMYD)
secMYD = CCI(pathL2SecMYD)

# ENVISAT AATSR paths and data
print "Reading ENVISAT AATSR data"
pathL2PriENV = mainL2 + "ESACCI-L2-CLOUD-CLD-AATSR_CC4CL_Envisat_200807221844_fv2.0.primary.nc"
pathL2SecENV = mainL1 + "ATS_TOA_1PUUPA20080722_184428_000065272070_00313_33433_6676.nc"
priMYD = CCI(pathL2PriMYD)
secMYD = CCI(pathL2SecMYD)

# N18 paths and data, RESAMPLED
print "Reading resampled N18 data"
pathL2N18PrimaryResampled = mainL2 + N18PrimaryResampledName
N18PrimaryResampled = CCI(pathL2N18PrimaryResampled)
pathL2N18SecondaryResampled = mainL2 + N18SecondaryResampledName
N18SecondaryResampled = CCI(pathL2N18SecondaryResampled)

# MODIS AQUA paths and data, RESAMPLED
print "Reading resampled MODIS AQUA data"
pathL2MYDPrimaryResampled = mainL2 + MYDPrimaryResampledName
MYDPrimaryResampled = CCI(pathL2MYDPrimaryResampled)
pathL2MYDSecondaryResampled = mainL2 + MYDSecondaryResampledName
MYDSecondaryResampled = CCI(pathL2MYDSecondaryResampled)

# ENVISAT AATSR paths and data, RESAMPLED
print "Reading resampled ENVISAT AATSR data"
pathL2ENVPrimaryResampled = mainL2 + ENVPrimaryResampledName
ENVPrimaryResampled = CCI(pathL2ENVPrimaryResampled)
pathL2ENVSecondaryResampled = mainL2 + ENVSecondaryResampledName
ENVSecondaryResampled = CCI(pathL2ENVSecondaryResampled)

# subset borders in lat/lon
centrePoint = [66.5, -108.]
leftPoint = [70, -150]
upperPoint = [87,  70]
rightPoint = [87.5,  82]
lowerPoint = [58, -85]
latBounds = [58, 87.5]
lonBounds = [-150,  95]
boundingBox = [leftPoint[1], upperPoint[1], lowerPoint[0], rightPoint[0]]
boundingBox[3] = 65.
boundingBox = [-134, -76, 52, 90]

# get all variables
MYDSlice = priMYD.getAllVariables(doSlice=True, boundingBox = boundingBox)
secMYD.getAllVariables(doSlice = True, boundingBox = boundingBox, primary = False, boxSlice = MYDSlice)
print "Getting all variables: N18 resampled"
N18PrimaryResampled.getAllVariables()
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
AllSensorsMaskCombined = N18ResampledCloudMask + MYDResampledCloudMask

# build mask of all pixels out of study area, i.e. where any sensor has no reflectance data
N18ReflMask = N18SecondaryResampled.reflectance_in_channel_no_1.mask
MYDReflMask = MYDSecondaryResampled.reflectance_in_channel_no_1.mask
ENVReflMask = ENVSecondaryResampled.reflectance_in_channel_no_1.mask
ReflMask = N18ReflMask + MYDReflMask + ENVReflMask

#################################################################################################
# plot variables and calculate statistics
# 1) compare CCI products from N18, MYD, and AATSR

    # a) RGB, ideally highly resolved MODIS with 0.6, 0.8, and 1.6 (1.6 has missing scan lines though)
if doRGB:    
    print "RGB: started."
    platform = "MYD"
    colourTupleMYD = buildRGB(MYDPrimaryResampled, MYDSecondaryResampled, platform)
    RGBName = figuresDir + "RGB_" + platform + "_" + delLatStr + "x" + delLonStr + ".png"
    plotRGB(RGBName, colourTupleMYD, MYDSecondaryResampled.lat, MYDSecondaryResampled.lon, MYDSecondaryResampled.reflectance_in_channel_no_1, 
            centrePoint[0], centrePoint[1])
    platform = "N18"
    colourTupleN18 = buildRGB(N18PrimaryResampled, N18SecondaryResampled, platform)
    RGBName = figuresDir + "RGB_" + platform + "_" + delLatStr + "x" + delLonStr + ".png"
    plotRGB(RGBName, colourTupleN18, N18SecondaryResampled.lat, N18SecondaryResampled.lon, N18SecondaryResampled.reflectance_in_channel_no_1, 
            centrePoint[0], centrePoint[1])
    platform = "ENV"
    colourTupleENV = buildRGB(ENVPrimaryResampled, ENVSecondaryResampled, platform)
    RGBName = figuresDir + "RGB_" + platform + "_" + delLatStr + "x" + delLonStr + ".png"
    plotRGB(RGBName, colourTupleENV, ENVSecondaryResampled.lat, ENVSecondaryResampled.lon, ENVSecondaryResampled.reflectance_in_channel_no_1, 
            centrePoint[0], centrePoint[1])    
    colourTupleMulti = np.concatenate((colourTupleN18[..., np.newaxis], colourTupleMYD[..., np.newaxis], colourTupleENV[..., np.newaxis]), axis=2)
    RGBName = figuresDir + "RGB_multi_" + delLatStr + "x" + delLonStr + ".png"
    plotRGBMulti(RGBName, colourTupleMulti, N18SecondaryResampled.lat, N18SecondaryResampled.lon, N18SecondaryResampled.reflectance_in_channel_no_1, 
                calipsoLat, calipsoLon,
                centrePoint[0], centrePoint[1])    
    print "RGB: done."    

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
collocateN18 = collocateCciAndCalipso(N18PrimaryResampled, calipsoData, maxDistance)
collocateMYD = collocateCciAndCalipso(MYDPrimaryResampled, calipsoData, maxDistance)
collocateENV = collocateCciAndCalipso(ENVPrimaryResampled, calipsoData, maxDistance)

"""Plot collocated data for COT, CTP, and CTT."""
plotCciCalipsoCollocation(collocateN18, collocateMYD, collocateENV, figuresDir)

N18PrimaryResampled.maskAllVariables(ReflMask)
N18SecondaryResampled.maskAllVariables(ReflMask)
MYDPrimaryResampled.maskAllVariables(ReflMask)
MYDSecondaryResampled.maskAllVariables(ReflMask)
ENVPrimaryResampled.maskAllVariables(ReflMask)
ENVSecondaryResampled.maskAllVariables(ReflMask)

sys.exit()

    
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

    # d) plot ctp for N18, MYD, AATSR

    # e) plot histogram with normal (?) distribution fit

    # f) calculate statistical moments: mean, median, std, skew, kurtosis
         # plot boxplots?

    # g) calculate significance of differences in variance and mean

    # h) calculate residuals and analyse dependence on third variable
        # RES = CTP_MYD - CTP_N18 -> mean/median/std
        # linear regression of RES with COT or TIMEDIFF
            # scatter plot RES vs COT/TIMEDIFF?

sys.exit()

print "Plotting data."
variable = "phase"

# 1a) RGB, ideally highly resolved MODIS with 0.6, 0.8, and 1.6 (1.6 has missing scan lines though) 

# 1b) calculate min/mean/max time difference between grid boxes = TIMEDIFF (no plot)

plotCCI(N18PrimaryResampled, MYDPrimaryResampled, boundingBox, centrePoint, variable,
        'NOAA18', mask = AllSensorsMaskCombined) 
print "N18 " + variable + " mean = " + str(round(getattr(N18PrimaryResampled, variable).mean(), 2))
print "N18 " + variable + " stdev = " + str(round(getattr(N18PrimaryResampled, variable).std(), 2))
plotCCI(N18PrimaryResampled, MYDPrimaryResampled, boundingBox, centrePoint, variable,
        'MYD', mask = AllSensorsMaskCombined) 
print "MYD " + variable + " mean = " + str(round(getattr(MYDPrimaryResampled, variable).mean(), 2))
print "MYD " + variable + " stdev = " + str(round(getattr(MYDPrimaryResampled, variable).std(), 2))
input = abs(getattr(MYDPrimaryResampled, variable) - getattr(N18PrimaryResampled, variable))
print [round(x, 2) for x in minMax(input)]
plotCCI(N18PrimaryResampled, MYDPrimaryResampled, boundingBox, centrePoint, variable,
        'MYDResampledMinusNOAA18', input = input, mask = input.mask) 
input.mask = input.mask + ma.masked_greater(input, input.mean() + input.std() * 2).mask
plotCCI(N18PrimaryResampled, MYDPrimaryResampled, boundingBox, centrePoint, variable,
        'MYDResampledMinusNOAA18', input = input, mask = input.mask) 
print "MYD-N18 " + variable + " mean = " + str(round(input.mean(), 2))
print "MYD-N18 " + variable + " stdev = " + str(round(input.std(), 2))
# plt.scatter(MYDPrimaryResampled.ctp, input)
# plt.draw()
# plt.scatter(timeDiff, input)
# plt.draw()
plt.show()

#stats.ttest_ind(getattr(N18PrimaryResampled, variable).ravel(), getattr(MYDPrimaryResampled, variable).ravel())
