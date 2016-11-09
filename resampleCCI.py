#!/data/osus/Enthought/User/bin/python2.7
# /home/oliver/Enthought/Canopy_64bit/User/bin/python

from analyseCCI import CCI, cciGrid
from CCITools import resampleCCI, writeCCI
import sys
from sys import argv

if len(argv) > 1:
    delLon = float(argv[1])
    delLat = float(argv[2])
    if argv[3] == "True":
        primary = True
    elif argv[3] == "False":
        primary = False
    else:
        print "ERROR: 3rd argument should be [True/False]."
        sys.exit()
else:
    delLat = 0.1
    delLon = 0.1
    primary = True

primaryString = "secondary"
if primary:
        primaryString = "primary"
print "Resampling " + primaryString + " data to " + str(delLon) + " lon x " + str(delLat) + " lat regular grid."

# main path to input files
mainL1 = "/cmsaf/cmsaf-cld7/esa_cloud_cci/data/v2.0/L1/"
mainL2 = "/cmsaf/cmsaf-cld7/esa_cloud_cci/data/v2.0/L2/"
if primary:
    suffix = "_primary"
else:
    suffix = "_secondary"
outNameN18 = "N18_Resampled_2008-07-22-1851_"
outNameMYD = "MYD_Resampled_2008-07-22-1915_"
outNameENV = "ENV_Resampled_2008-07-22-1844_"

# targetGrid:
minLat = 40.0
maxLat = 90.0
minLon = -179.0
maxLon = 0.

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
boundingBox = [-179., 0., 40, 90] #[-134, -76, 52, 85]

# NOAA18 paths and data
print "Reading NOAA18 data"
pathL2PriN18 = mainL2 + "20080722185100-ESACCI-L2_CLOUD-CLD_PRODUCTS-AVHRRGAC-NOAA18-fv2.0.nc"
pathL2SecN18 = mainL2 + "ECC_GAC_avhrr_noaa18_99999_20080722T1851289Z_20080722T2046134Z.secondary.nc"
priN18 = CCI(pathL2PriN18)
secN18 = CCI(pathL2SecN18)

# MODIS AQUA paths and data
print "Reading MODIS AQUA data"
pathL2PriMYD = mainL2 + "MYD_merged_20080722_19151920_primary.nc" #"MYD20080722_1915.nc"
pathL2SecMYD = mainL2 + "MYD_merged_20080722_19151920_secondary.nc" # "MYD021KM.A2008204.1915.006.2012069115248.bspscs_000500694537.secondary.nc"
priMYD = CCI(pathL2PriMYD)
secMYD = CCI(pathL2SecMYD)
# pathL2PriMYD2 = mainL2 + "MYD20080722_1920.nc"
# pathL2SecMYD2 = mainL2 + "MYD021KM.A2008204.1920.006.2012069113627.bspscs_000500694537.secondary.nc"
# mergeGranules(pathL2PriMYD, pathL2PriMYD2, mainL2 + "MYD_merged_20080722_19151920_primary.nc")
# mergeGranules(pathL2SecMYD, pathL2SecMYD2, mainL2 + "MYD_merged_20080722_19151920_secondary.nc")
# priMYD2 = CCI(pathL2PriMYD2)
# secMYD2 = CCI(pathL2SecMYD2)

# Envisat AATSR paths and data
print "Reading ENVISAT AATSR data"
pathL2PriENV = mainL2 + "ESACCI-L2-CLOUD-CLD-AATSR_CC4CL_Envisat_200807221844_fv2.0.primary.nc"
pathL2SecENV = mainL1 + "ATS_TOA_1PUUPA20080722_184428_000065272070_00313_33433_6676.nc"
priENV = CCI(pathL2PriENV)
secENV = CCI(pathL2SecENV)

# get all variables
print "Getting all variables: NOAA18"
N18Slice = priN18.getAllVariables(doSlice = True, boundingBox = boundingBox)
secN18.getAllVariables(doSlice = True, boundingBox = boundingBox, primary = False, boxSlice = N18Slice)
print "Getting all variables: MODIS AQUA"
MYDSlice = priMYD.getAllVariables(doSlice = True, boundingBox = boundingBox)
secMYD.getAllVariables(doSlice = True, boundingBox = boundingBox, primary = False, boxSlice = MYDSlice)
# MYDSlice2 = priMYD2.getAllVariables(doSlice = True, boundingBox = boundingBox)
# secMYD2.getAllVariables(doSlice = True, boundingBox = boundingBox, primary = False, boxSlice = MYDSlice2)
print "Getting all variables: ENVISAT AATSR"
ENVSlice = priENV.getAllVariables(doSlice = True, boundingBox = boundingBox)
ENVSlice = secENV.getAllVariables(doSlice = True, boundingBox = boundingBox)
secENV.reflec_nadir_0670 /= secENV.reflec_nadir_0670.max()
secENV.reflec_nadir_0870 /= secENV.reflec_nadir_0870.max()
secENV.reflec_nadir_1600 /= secENV.reflec_nadir_1600.max()

# resample to N18 if requested
print "Resampling data."
# primary
if primary:
    resampledData = resampleCCI(priMYD, targetGrid, "MYD")
    fileName = mainL2 + outNameMYD + "lat" + str(delLat) + "lon" + str(delLon) + suffix + ".nc"
    print "    writing file " + fileName
    writeCCI(fileName, resampledData, targetGrid, primary)
    resampledData = resampleCCI(priN18, targetGrid, "N18")
    fileName = mainL2 + outNameN18 + "lat" + str(delLat) + "lon" + str(delLon) + suffix + ".nc"
    print "    writing file " + fileName
    writeCCI(fileName, resampledData, targetGrid, primary)
    resampledData = resampleCCI(priENV, targetGrid, "ENV")
    fileName = mainL2 + outNameENV + "lat" + str(delLat) + "lon" + str(delLon) + suffix + ".nc"
    print "    writing file " + fileName
    writeCCI(fileName, resampledData, targetGrid, primary)
# secondary
else:
    resampledData = resampleCCI(secMYD, targetGrid, "MYD", lat_in = priMYD.lat, lon_in = priMYD.lon) 
    fileName = mainL2 + outNameMYD + "lat" + str(delLat) + "lon" + str(delLon) + suffix + ".nc"
    print "    writing file " + fileName
    writeCCI(fileName, resampledData, targetGrid, primary, platform = "MYD")
    resampledData = resampleCCI(secN18, targetGrid, "N18", lat_in = priN18.lat, lon_in = priN18.lon) 
    fileName = mainL2 + outNameN18 + "lat" + str(delLat) + "lon" + str(delLon) + suffix + ".nc"
    print "    writing file " + fileName
    writeCCI(fileName, resampledData, targetGrid, primary, platform = "N18")
    resampledData = resampleCCI(secENV, targetGrid, "ENV", lat_in = secENV.lat, lon_in = secENV.lon)
    fileName = mainL2 + outNameENV + "lat" + str(delLat) + "lon" + str(delLon) + suffix + ".nc"
    print "    writing file " + fileName
    writeCCI(fileName, resampledData, targetGrid, primary, platform = "ENV")






