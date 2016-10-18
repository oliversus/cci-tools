#!/cmsaf/nfshome/routcm/Modules_sw/python/2.7.9/bin/python

from analyseCCI import CCI
from mpl_toolkits.basemap import Basemap, cm
import netCDF4
import matplotlib.pyplot as plt
from matplotlib import path
import numpy as np
from CCITools import minMax, plotCCI, buildRGB, plotRGB
import sys
import math
import numpy.ma as ma
from scipy import stats
from osgeo._gdal import GDAL_GCP_GCPLine_get

figuresDir = "/cmsaf/esa_doku/ESA_Cloud_cci/publications/CC4CL_paper/figures/"
mainL2 = "/cmsaf/cmsaf-cld7/esa_cloud_cci/data/v2.0/L2/"
delLat = "0.05"
delLon = "0.05"
N18PrimaryResampledName = "N18_Resampled_2008-07-22-1851_lat" + delLat + "lon" + delLon + "_primary.nc"
MYDPrimaryResampledName = "MYD_Resampled_2008-07-22-1915_lat" + delLat + "lon" + delLon + "_primary.nc"
N18SecondaryResampledName = "N18_Resampled_2008-07-22-1851_lat" + delLat + "lon" + delLon + "_secondary.nc"
MYDSecondaryResampledName = "MYD_Resampled_2008-07-22-1915_lat" + delLat + "lon" + delLon + "_secondary.nc"

# NOAA18 paths and data
print "Reading NOAA18 data"
pathL2PriN18 = mainL2 + "20080722185100-ESACCI-L2_CLOUD-CLD_PRODUCTS-AVHRRGAC-NOAA18-fv2.0.nc"
pathL2SecN18 = mainL2 + "ECC_GAC_avhrr_noaa18_99999_20080722T1851289Z_20080722T2046134Z.secondary.nc"   
priN18 = CCI(pathL2PriN18)
secN18 = CCI(pathL2SecN18)

# MODIS AQUA paths and data
print "Reading MODIS AQUA data"
pathL2PriMYD = mainL2 + "MYD20080722_1915.nc"
pathL2SecMYD = mainL2 + "MYD021KM.A2008204.1915.006.2012069115248.bspscs_000500694537.secondary.nc"
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

#ds = gdal.Open(pathL2MYDSecondaryResampled).ReadAsArray()


# subset borders in lat/lon
centrePoint = [64.5, -102.5] #[71 ,  -80 ] 
leftPoint   = [70  , -150]
upperPoint  = [87  ,  70 ]#95 ]
rightPoint  = [87.5,  82 ]
lowerPoint  = [58  , -85 ]
latBounds   = [58  , 87.5]
lonBounds   = [-150,  95 ]
boundingBox = [leftPoint[1], upperPoint[1], lowerPoint[0], rightPoint[0]]
boundingBox[3] = 65.
boundingBox = [-134, -76, 52, 74]

# get all variables
# print "Getting all variables: NOAA18"
# N18Slice = priN18.getAllVariables(doSlice = True, boundingBox = boundingBox)
# secN18.getAllVariables(doSlice = True, boundingBox = boundingBox, primary = False, boxSlice = N18Slice)
# print "Getting all variables: MODIS AQUA"
MYDSlice = priMYD.getAllVariables(doSlice = True, boundingBox = boundingBox)
secMYD.getAllVariables(doSlice = True, boundingBox = boundingBox, primary = False, boxSlice = MYDSlice)
print "Getting all variables: N18 resampled"
N18PrimaryResampled.getAllVariables()
N18SecondaryResampled.getAllVariables()
print "Getting all variables: MODIS resampled"
MYDPrimaryResampled.getAllVariables()
MYDSecondaryResampled.getAllVariables()

# mask all resampled pixels with cc_total < 1 to exclude fractional cloud coverage
N18ResampledCloudMask = ma.masked_less(N18PrimaryResampled.cc_total, 1.).mask
MYDResampledCloudMask = ma.masked_less(MYDPrimaryResampled.cc_total, 1.).mask
N18MYDMaskCombined = N18ResampledCloudMask + MYDResampledCloudMask

#################################################################################################
# plot variables and calculate statistics
# 1) compare CCI products from N18, MYD, and AATSR

    # a) RGB, ideally highly resolved MODIS with 0.6, 0.8, and 1.6 (1.6 has missing scan lines though)
if False:
    print "RGB: started."
    colourTuple = buildRGB(MYDPrimaryResampled, MYDSecondaryResampled)
    plotRGB(colourTuple, MYDSecondaryResampled.lat, MYDSecondaryResampled.lon, MYDSecondaryResampled.albedo_in_channel_no_1, 
            centrePoint[0], centrePoint[1])
    print "RGB: done."

    # b) calculate min/mean/max time difference between grid boxes = TIMEDIFF (no plot)

    # c) calculate min/mean/max solzen and satzen for all 3 sensors

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

print "Plotting data."
variable = "cth"

# 1a) RGB, ideally highly resolved MODIS with 0.6, 0.8, and 1.6 (1.6 has missing scan lines though) 

# 1b) calculate min/mean/max time difference between grid boxes = TIMEDIFF (no plot)

timeDiff = ((MYDResampled.time - math.floor(MYDResampled.time.min())) \
                - (N18Resampled.time - math.floor(N18Resampled.time.min()))) * 24 * 60
if True:
    print "The minimum/maximum time difference between MODIS AQUA and NOAA18 is: ", \
        round(timeDiff.min(), 2), "min /", round(timeDiff.max(), 2), "min"



plotCCI(N18Resampled, MYDResampled, boundingBox, centrePoint, variable, 
        'NOAA18', mask = N18MYDMaskCombined) 
print "N18 " + variable + " mean = " + str(round(getattr(N18Resampled, variable).mean(), 2))
print "N18 " + variable + " stdev = " + str(round(getattr(N18Resampled, variable).std(), 2))
plotCCI(N18Resampled, MYDResampled, boundingBox, centrePoint, variable, 
        'MYD', mask = N18MYDMaskCombined) 
print "MYD " + variable + " mean = " + str(round(getattr(MYDResampled, variable).mean(), 2))
print "MYD " + variable + " stdev = " + str(round(getattr(MYDResampled, variable).std(), 2))
input = abs(getattr(MYDResampled, variable) - getattr(N18Resampled, variable))
print [round(x, 2) for x in minMax(input)]
plotCCI(N18Resampled, MYDResampled, boundingBox, centrePoint, variable, 
        'MYDResampledMinusNOAA18', input = input, mask = input.mask) 
input.mask = input.mask + ma.masked_greater(input, input.mean() + input.std() * 2).mask
plotCCI(N18Resampled, MYDResampled, boundingBox, centrePoint, variable, 
        'MYDResampledMinusNOAA18', input = input, mask = input.mask) 
print "MYD-N18 " + variable + " mean = " + str(round(input.mean(), 2))
print "MYD-N18 " + variable + " stdev = " + str(round(input.std(), 2))
# plt.scatter(MYDResampled.ctp, input)
# plt.draw()
# plt.scatter(timeDiff, input)
# plt.draw()
plt.show()

stats.ttest_ind(getattr(N18Resampled, variable).ravel(), getattr(MYDResampled, variable).ravel())
