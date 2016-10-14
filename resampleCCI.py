#!/cmsaf/nfshome/routcm/Modules_sw/python/2.7.9/bin/python

from analyseCCI import CCI, cciGrid
from mpl_toolkits.basemap import Basemap, cm
import netCDF4
import matplotlib.pyplot as plt
from matplotlib import path
import numpy as np
from CCITools import resampleCCI, minMax, writeCCI, plotCCI
import sys
import math
import numpy.ma as ma

# main path to input files
mainL2 = "/cmsaf/cmsaf-cld7/esa_cloud_cci/data/v2.0/L2/"
primary = True
if primary:
    suffix = "_primary"
else:
    suffix = "_secondary"
outNameN18 = "N18_Resampled_2008-07-22-1851_"
outNameMYD = "MYD_Resampled_2008-07-22-1915_"

# targetGrid:
minLat = 52.0
maxLat = 74.0
delLat = 0.05
minLon = -137.0
maxLon = -76.0
delLon = 0.05

targetGrid = cciGrid(minLat, minLon, maxLat, maxLon, delLat = delLat, delLon = delLon)
targetGrid.buildGrid()

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

# get all variables
print "Getting all variables: NOAA18"
N18Slice = priN18.getAllVariables(doSlice = True, boundingBox = boundingBox)
secN18.getAllVariables(doSlice = True, boundingBox = boundingBox, primary = False, boxSlice = N18Slice)
print "Getting all variables: MODIS AQUA"
MYDSlice = priMYD.getAllVariables(doSlice = True, boundingBox = boundingBox)
secMYD.getAllVariables(doSlice = True, boundingBox = boundingBox, primary = False, boxSlice = MYDSlice)

# resample to N18 if requested
print "Resampling data."
# primary
if primary:
    resampledData = resampleCCI(priMYD, targetGrid)
    fileName = mainL2 + outNameMYD + "lat" + str(delLat) + "lon" + str(delLon) + suffix + ".nc"
    print "    writing file " + fileName
    writeCCI(fileName, resampledData, targetGrid, primary)
    resampledData = resampleCCI(priN18, targetGrid)
    fileName = mainL2 + outNameN18 + "lat" + str(delLat) + "lon" + str(delLon) + suffix + ".nc"
    print "    writing file " + fileName
    writeCCI(fileName, resampledData, targetGrid, primary)
# secondary
else:
    resampledData = resampleCCI(secMYD, targetGrid, lat_in = priMYD.lat, lon_in = priMYD.lon, NOAA = False) 
    fileName = mainL2 + outNameMYD + "lat" + str(delLat) + "lon" + str(delLon) + suffix + ".nc"
    print "    writing file " + fileName
    writeCCI(fileName, resampledData, targetGrid, primary, NOAA = False)
    resampledData = resampleCCI(secN18, targetGrid, lat_in = priN18.lat, lon_in = priN18.lon, NOAA = True) 
    fileName = mainL2 + outNameN18 + "lat" + str(delLat) + "lon" + str(delLon) + suffix + ".nc"
    print "    writing file " + fileName
    writeCCI(fileName, resampledData, targetGrid, primary, NOAA = True)
