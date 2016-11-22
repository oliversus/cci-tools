#!/data/osus/Enthought/User/bin/python2.7

from analyseCCI import CCI, cciGrid
from CCITools import resampleCCI, writeCCI, greatCircle
import sys
from sys import argv
import os

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
    sceneTime = str(argv[4])  # 1) 07221915 2) 07270810
    if sceneTime != '07221915' and sceneTime != '07270810' and sceneTime != '07230021' and sceneTime != '07222058':
        print "ERROR: choose correct study date ('07221915' or '07270810' or '07230021' or '07222058')"
        sys.exit()
else:
    delLat = 0.1
    delLon = 0.1
    primary = True
    sceneTime = '07222058'

month = sceneTime[0:2]
day = sceneTime[2:4]
hour = sceneTime[4:6]
minute = sceneTime[6:8]

data_folder = "/cmsaf/esa_doku/ESA_Cloud_cci/publications/CC4CL_paper/data/"

primaryString = "secondary"
if primary:
        primaryString = "primary"
print "Resampling " + primaryString + " data to " + str(delLon) + " lon x " + str(delLat) + " lat regular grid."

if primary:
    suffix = "_primary"
else:
    suffix = "_secondary"

outNameN18 = "N18_Resampled_2008-" + month + "-" + day + "-" + hour + minute + "_"
outNameMYD = outNameN18.replace("N18", "MYD")
outNameENV = outNameN18.replace("N18", "ENV") #"ENV" + outNameN18[3:len(outNameN18)]

# subset borders in lat/lon
if sceneTime == '07221915':
    centrePoint = [64.5, -102.5]
    boundingBox = [-179., 0., 40, 90]
elif sceneTime == '07270810':
    centrePoint = [73., 55.]
    boundingBox = [-10., 130., 45., 90.]
elif sceneTime == '07230021':
    centrePoint = [71., 173.]
    boundingBox = [140., -160., 45., 90.]
elif sceneTime == '07222058':
    centrePoint = [74., -143.]
    boundingBox = [-180., -100., 45., 90.]

# targetGrid:
minLat = boundingBox[2]
maxLat = boundingBox[3]
minLon = boundingBox[0]
maxLon = boundingBox[1]

targetGrid = cciGrid(minLat, minLon, maxLat, maxLon, delLat = delLat, delLon = delLon)
targetGrid.buildGrid()

# NOAA18 paths and data
l2_primary_prefix = "cci_l2_primary_"
print "Reading NOAA18 data"
pathL2PriN18 = data_folder + l2_primary_prefix + "n18_" + sceneTime + ".nc"
pathL2SecN18 = pathL2PriN18.replace("primary", "secondary")
priN18 = CCI(pathL2PriN18)
secN18 = CCI(pathL2SecN18)

# MODIS AQUA paths and data
print "Reading MODIS AQUA data"
pathL2PriMYD = data_folder + l2_primary_prefix + "myd_" + sceneTime + ".nc"
pathL2SecMYD = pathL2PriMYD.replace("primary", "secondary")
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
pathL2PriENV = data_folder + l2_primary_prefix + "env_" + sceneTime + ".nc"
pathL2SecENV = pathL2PriENV.replace("primary", "secondary")
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

"""the maximum distance between CCI and a grid box center is given by
    the great circle distance between the grid box center and one of its corners
    which is half of the grid box diagonal, which can be calculated with Pythagoras"""
"""The lat/lon spacing and the grid box longitude/latitude width"""
boxLonWidth = delLon * greatCircle(centrePoint[1], centrePoint[0], centrePoint[1] + 1., centrePoint[0])
boxLatWidth = delLat * greatCircle(centrePoint[1], centrePoint[0], centrePoint[1], centrePoint[0] + 1.)
"""give the grid box diagonal"""
boxDiag = (boxLonWidth**2 + boxLatWidth**2)**0.5
"""half of which is the maximum distance of CCI to any grid box center"""
maxDistance = boxDiag / 2.
"""add 10 percent to this distance so that """
maxDistance *= 1.1

# resample to N18 if requested
print "Resampling data."
if primary:
    resampledData = resampleCCI(priMYD, targetGrid, "MYD", maxDistance)
    fileName = os.path.splitext(pathL2PriMYD)[0] + "_resampled_" + str(delLon) + "_" + str(delLat) + ".nc"
    print "    writing file " + fileName
    writeCCI(fileName, resampledData, targetGrid, primary)
    resampledData = resampleCCI(priN18, targetGrid, "N18", maxDistance)
    fileName = os.path.splitext(pathL2PriN18)[0] + "_resampled_" \
                                                   "" \
                                                   "" + str(delLon) + "_" + str(delLat) + ".nc"
    print "    writing file " + fileName
    writeCCI(fileName, resampledData, targetGrid, primary)
    resampledData = resampleCCI(priENV, targetGrid, "ENV", maxDistance)
    fileName = os.path.splitext(pathL2PriENV)[0] + "_resampled_" + str(delLon) + "_" + str(delLat) + ".nc"
    print "    writing file " + fileName
    writeCCI(fileName, resampledData, targetGrid, primary)
else:
    resampledData = resampleCCI(secMYD, targetGrid, "MYD", maxDistance, lat_in = priMYD.lat, lon_in = priMYD.lon)
    fileName = os.path.splitext(pathL2SecMYD)[0] + "_resampled_" + str(delLon) + "_" + str(delLat) + ".nc"
    print "    writing file " + fileName
    writeCCI(fileName, resampledData, targetGrid, primary, platform = "MYD")
    resampledData = resampleCCI(secN18, targetGrid, "N18", maxDistance, lat_in = priN18.lat, lon_in = priN18.lon)
    fileName = os.path.splitext(pathL2SecN18)[0] + "_resampled_" + str(delLon) + "_" + str(delLat) + ".nc"
    print "    writing file " + fileName
    writeCCI(fileName, resampledData, targetGrid, primary, platform = "N18")
    resampledData = resampleCCI(secENV, targetGrid, "ENV", maxDistance, lat_in = secENV.lat, lon_in = secENV.lon)
    fileName = os.path.splitext(pathL2SecENV)[0] + "_resampled_" + str(delLon) + "_" + str(delLat) + ".nc"
    print "    writing file " + fileName
    writeCCI(fileName, resampledData, targetGrid, primary, platform = "ENV")






