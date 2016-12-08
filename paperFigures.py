#!/data/osus/Enthought/User/bin/python2.7

from analyseCCI import CCI
import numpy as np
from CCITools import buildRGB, plotRGB,\
    plotRGBMulti, greatCircle, collocateCciAndCalipso, \
    plotCciCalipsoCollocation, plotCCI, update_latex_variables
import sys
import numpy.ma as ma
from pyhdf.SD import SD, SDC
from sys import argv
import math
import globals

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
    if argv[4] != '07221915' and argv[4] != '07270810' and argv[4] != '07230021' and argv[4] != '07222058':
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
    delLat = "0.5"
    delLon = "0.5"
    doRGB = False
    sceneTime = '07221915'
    corrected = True
    plotCot = False
    plotCalipso = True

delLat_delLon = delLat + "_" + delLon
month = sceneTime[0:2]
day = sceneTime[2:4]
hour = sceneTime[4:6]
minute = sceneTime[6:8]

data_folder = "/cmsaf/esa_doku/ESA_Cloud_cci/publications/CC4CL_paper/data/"
figuresDir = "/cmsaf/esa_doku/ESA_Cloud_cci/publications/CC4CL_paper/figures/"

calipsoPath1km = data_folder + "calipso_1km_" + sceneTime + ".hdf"
calipsoPath5km = data_folder + "calipso_5km_" + sceneTime + ".hdf"

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

# ENVISAT AATSR paths and data
print "Reading ENVISAT AATSR data"
pathL2PriENV = data_folder + l2_primary_prefix + "env_" + sceneTime + ".nc"
pathL2SecENV = pathL2PriENV.replace("primary", "secondary")
priENV = CCI(pathL2PriENV)
secENV = CCI(pathL2SecENV)

# N18 paths and data, RESAMPLED
print "Reading resampled N18 data"
pathL2N18PrimaryResampled = data_folder + N18PrimaryResampledName
N18PrimaryResampled = CCI(pathL2N18PrimaryResampled)
N18Masked = CCI(pathL2N18PrimaryResampled)
pathL2N18SecondaryResampled = data_folder + N18SecondaryResampledName
N18SecondaryResampled = CCI(pathL2N18SecondaryResampled)

# MODIS AQUA paths and data, RESAMPLED
print "Reading resampled MODIS AQUA data"
pathL2MYDPrimaryResampled = data_folder + MYDPrimaryResampledName
MYDPrimaryResampled = CCI(pathL2MYDPrimaryResampled)
pathL2MYDSecondaryResampled = data_folder + MYDSecondaryResampledName
MYDSecondaryResampled = CCI(pathL2MYDSecondaryResampled)

# ENVISAT AATSR paths and data, RESAMPLED
print "Reading resampled ENVISAT AATSR data"
pathL2ENVPrimaryResampled = data_folder + ENVPrimaryResampledName
ENVPrimaryResampled = CCI(pathL2ENVPrimaryResampled)
pathL2ENVSecondaryResampled = data_folder + ENVSecondaryResampledName
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
MYDSlice = priMYD.getAllVariables(doSlice=True, boundingBox=boundingBox)
secMYD.getAllVariables(doSlice=True, boundingBox=boundingBox, primary=False, boxSlice=MYDSlice)
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
AllSensorsMaskCombined = N18ResampledCloudMask + MYDResampledCloudMask

# build mask of all pixels out of study area, i.e. where any sensor has no reflectance data
N18ReflMask = N18SecondaryResampled.reflectance_in_channel_no_1.mask
MYDReflMask = MYDSecondaryResampled.reflectance_in_channel_no_1.mask
ENVReflMask = ENVSecondaryResampled.reflectance_in_channel_no_1.mask
ReflMask = N18ReflMask + MYDReflMask + ENVReflMask

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
    RGBName = figuresDir + "RGB_" + platform + "_" + delLatStr + "x" + delLonStr + "_" + sceneTime + ".png"
    plotRGB(RGBName, colourTupleMYD, MYDSecondaryResampled.lat, MYDSecondaryResampled.lon, MYDSecondaryResampled.reflectance_in_channel_no_1,
            centrePoint[0], centrePoint[1])
    platform = "N18"
    colourTupleN18 = buildRGB(N18PrimaryResampled, N18SecondaryResampled, platform)
    RGBName = figuresDir + "RGB_" + platform + "_" + delLatStr + "x" + delLonStr + "_" + sceneTime + ".png"
    plotRGB(RGBName, colourTupleN18, N18SecondaryResampled.lat, N18SecondaryResampled.lon, N18SecondaryResampled.reflectance_in_channel_no_1,
            centrePoint[0], centrePoint[1])
    platform = "ENV"
    colourTupleENV = buildRGB(ENVPrimaryResampled, ENVSecondaryResampled, platform)
    RGBName = figuresDir + "RGB_" + platform + "_" + delLatStr + "x" + delLonStr + "_" + sceneTime + ".png"
    plotRGB(RGBName, colourTupleENV, ENVSecondaryResampled.lat, ENVSecondaryResampled.lon, ENVSecondaryResampled.reflectance_in_channel_no_1,
            centrePoint[0], centrePoint[1])
    colourTupleMulti = np.concatenate((colourTupleN18[..., np.newaxis], colourTupleMYD[..., np.newaxis], colourTupleENV[..., np.newaxis]), axis=2)
    RGBName = figuresDir + "RGB_multi_" + delLatStr + "x" + delLonStr + "_" + sceneTime + ".png"
    plotRGBMulti(RGBName, colourTupleMulti, N18SecondaryResampled.lat, N18SecondaryResampled.lon, poly_lats, poly_lons, N18SecondaryResampled.reflectance_in_channel_no_1,
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
collocateN18 = collocateCciAndCalipso(N18PrimaryResampled, calipsoData, maxDistance, corrected)
collocateMYD = collocateCciAndCalipso(MYDPrimaryResampled, calipsoData, maxDistance, corrected)
collocateENV = collocateCciAndCalipso(ENVPrimaryResampled, calipsoData, maxDistance, corrected)

"""Plot collocated data for COT, CTP, and CTT."""
figurePath = figuresDir + "calipsoVsCci_" + sceneTime
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

N18PrimaryResampled.maskAllVariables(ReflMask)
N18SecondaryResampled.maskAllVariables(ReflMask)
MYDPrimaryResampled.maskAllVariables(ReflMask)
MYDSecondaryResampled.maskAllVariables(ReflMask)
ENVPrimaryResampled.maskAllVariables(ReflMask)
ENVSecondaryResampled.maskAllVariables(ReflMask)

#sys.exit()
    
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

#sys.exit()

print "Plotting data."

# 1a) RGB, ideally highly resolved MODIS with 0.6, 0.8, and 1.6 (1.6 has missing scan lines though) 

# 1b) calculate min/mean/max time difference between grid boxes = TIMEDIFF (no plot)
variable = "lat"
plotCCI(N18PrimaryResampled, MYDPrimaryResampled, boundingBox, centrePoint, variable,
        'NOAA18', mask = ReflMask)
print "N18 " + variable + " mean = " + str(round(getattr(N18PrimaryResampled, variable).mean(), 2))
print "N18 " + variable + " stdev = " + str(round(getattr(N18PrimaryResampled, variable).std(), 2))

variable = "lon"
plotCCI(N18PrimaryResampled, MYDPrimaryResampled, boundingBox, centrePoint, variable,
        'NOAA18', mask = ReflMask)
print "N18 " + variable + " mean = " + str(round(getattr(N18PrimaryResampled, variable).mean(), 2))
print "N18 " + variable + " stdev = " + str(round(getattr(N18PrimaryResampled, variable).std(), 2))

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

#stats.ttest_ind(getattr(N18PrimaryResampled, variable).ravel(), getattr(MYDPrimaryResampled, variable).ravel())

print "updating latex variables"
path = globals.main_folder
update_latex_variables(path)
print "...done"

# """loop over all tex files in main folder"""
# for tex in glob.glob(path + "*.tex"):
#     """open a test version of each for writing"""
#     with open(os.path.splitext(tex)[0] + "_test.tex", 'w') as new_f:
#         """open tex file for reading"""
#         with open(tex, 'r') as f:
#             """loop over all lines of each file"""
#             for line in f:
#                 """and check whether that line matches the search pattern"""
#                 if re.search('insertVariable{.*}[0-9]*\.{1}[0-9]*', line):
#                     """if so, loop over all words in line"""
#                     for word in line.split():
#                         """and check whether that word matches the search pattern"""
#                         m = re.search('insertVariable{.*}[0-9]*\.{1}[0-9]*', word)
#                         if m:
#                             """if so, get the matching word"""
#                             found = m.group()
#                             """now loop over all latex variables in dictionary"""
#                             for variable, value in globals.latex_variables.iteritems():
#                                 print variable, value
#                                 print variable in found
#                                 """and if a key matches the word"""
#                                 if variable in found:
#                                     value_new = str(value)
#                                     """get the old value of the latex variable from the tex file"""
#                                     value_old = re.search('[0-9]*\.[0-9]*', found).group()
#                                     """replace it with the new value"""
#                                     replace = found.replace(value_old, value_new)
#                                     """and replace the old word with that new word within the entire line"""
#                                     line = line.replace(found, replace)
#                 """write each line to the test output file, regardless of whether a match has been found"""
#                 new_f.write(line)
#     """after having looped over all lines of a file, replace it by its test version"""
#     os.remove(tex)
#     os.rename(os.path.splitext(tex)[0] + "_test.tex", tex)