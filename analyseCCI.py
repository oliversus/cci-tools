import netCDF4
from matplotlib import path
import numpy as np

class CCI():

    def __init__(self, path):
        self.path = path
        self.dataset = netCDF4.Dataset(self.path)

    def __repr__(self):
        return "CCI: " + str(self.path)

    def openNetCDF(self):
        self.dataset = netCDF4.Dataset(self.path)

    def getPath(self):
        return self.path

    def setPath(self, path):
        self.path = path

    def printAllVariables(self):
        for i in iter(self.dataset.variables):
            print i + ",",

    def getVariable(self, vName, slice = None):
        if slice:
            setattr(self, vName, self.dataset.variables[vName][slice])
        else:
            setattr(self, vName, self.dataset.variables[vName][:])

    def getAllVariables(self, doSlice = False, boundingBox = [-150, 70, 58, 87.5], primary = True, boxSlice = []):
        if doSlice:
            if primary:
                self.getSlice(boundingBox = boundingBox)
                for vName in iter(self.dataset.variables):
                    setattr(self, vName, self.dataset.variables[vName][self.boxSlice])
                return self.boxSlice
            else:
                self.boxSlice = boxSlice
                for vName in iter(self.dataset.variables):
                    setattr(self, vName, self.dataset.variables[vName][self.boxSlice])
        else:
            for vName in iter(self.dataset.variables):
                setattr(self, vName, self.dataset.variables[vName][:])

    def getLon(self, slice = None):
        if slice:
            self.lon = self.dataset.variables['lon'][slice]
        else:
            self.lon = self.dataset.variables['lon'][:]
        return

    def getLat(self, slice = None):
        if slice:
            self.lat = self.dataset.variables['lat'][slice]
        else:
            self.lat = self.dataset.variables['lat'][:]
        return

    def getLatLon(self, slice = None):
        if slice:
            self.lat = self.dataset.variables['lat'][slice]
            self.lon = self.dataset.variables['lon'][slice]
        else:
            self.lat = self.dataset.variables['lat'][:]
            self.lon = self.dataset.variables['lon'][:]
        return

    def minmax(self, a):
        return a.min(), a.max()

    def getBoundingBoxIndexes(self, bbox = [-150, 70, 58, 87.5]):
        try:
            bbox = np.array(bbox)
            mypath = np.array([bbox[[0, 1, 1, 0]], bbox[[2, 2, 3, 3]]]).T
            p = path.Path(mypath)
            points = np.vstack((self.lon.flatten(), self.lat.flatten())).T   
            n, m = np.shape(self.lon)
            inside = p.contains_points(points).reshape((n, m))
            ii, jj = np.meshgrid(xrange(m), xrange(n))
            return [min(ii[inside]), max(ii[inside]), min(jj[inside]), max(jj[inside])]
        except AttributeError as detail: 
            print "Attribute error: ", detail

    def getLatLonSliced(self, boundingBox = [-150, 70, 58, 65.0]):
        self.getLatLon()
        boundingBoxIndexes = self.getBoundingBoxIndexes(boundingBox)
        self.boxSlice = np.index_exp[boundingBoxIndexes[2]:boundingBoxIndexes[3],
                                     boundingBoxIndexes[0]:boundingBoxIndexes[1]]
        self.getLatLon(slice = self.boxSlice)

    def getSlice(self, boundingBox = [-150, 70, 58, 65.0]):
        self.getLatLon()
        boundingBoxIndexes = self.getBoundingBoxIndexes(boundingBox)
        self.boxSlice = np.index_exp[boundingBoxIndexes[2]:boundingBoxIndexes[3],
                                     boundingBoxIndexes[0]:boundingBoxIndexes[1]]
        del self.lon
        del self.lat
        
class cciGrid():
    
    def __init__(self, minLat, minLon, maxLat, maxLon, delLat = 0.1, delLon = 0.1):
        self.minLat = minLat
        self.minLon = minLon
        self.maxLat = maxLat
        self.maxLon = maxLon
        self.delLat = delLat
        self.delLon = delLon        

    def __repr__(self):
        return "CCI grid min/max values: lat = " + str(self.minLat) + "/" + str(self.maxLat)+ \
            ", lon = " + str(self.minLon) + "/" + str(self.maxLon) 

    def buildGrid(self):
        self.latIter = self.buildVector(self.minLat, self.maxLat, self.delLat)
        self.lonIter = self.buildVector(self.minLon, self.maxLon, self.delLon)

        self.latVector = np.fromiter(self.latIter, dtype = np.float)
        self.lonVector = np.fromiter(self.lonIter, dtype = np.float)

        self.lat, self.lon = np.meshgrid(self.latVector, self.lonVector ,indexing = 'xy')

    def buildVector(self, start, end, stepSize):
        val = start
        if end < start:
            end += 360.0
        while val < end + stepSize:
            if (val > 180.0):
                yield val - 360.0
            else:
                yield val
            val += stepSize
        
        

    

