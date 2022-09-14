# -*- coding: utf-8 -*-
"""
Created on Aug 1 00:00:00 2022
Based on the code written by jlwoodr3

@author: shinbunya
"""

#from turtle import end_fill
#from zmq import curve_keypair
import numpy as np
#import cupy as cp


class Control:
    """Class to store info in control file"""
    def __init__(self):
        self.loaded = False

class Mesh:
    """Class to store info in mesh file"""
    def __init__(self):
        self.loaded = False

class Subgrid:
    """Class to store subgrid data"""
    def __init__(self):
        self.loaded = False

class SubgridVectorized:
    """Class to store vectorized subgrid data"""
    def __init__(self):
        self.loaded = False

class SubgridDepthBased:
    """Class to store depth-based subgrid data"""
    def __init__(self):
        self.loaded = False

class SubgridCalculatorDGP0():
    """Class to evaluate and hold subgrid data"""

    def __init__(self, subgridControlFilename, level0andLevel1 = True, GPU = True, \
        numHBasedTableLevels = 11):
        import numpy

        self.level0andLevel1 = level0andLevel1  # state if you want level0 and level1 corrections or just level 0
        self.control = Control()
        self.mesh = Mesh()
        self.subgrid = Subgrid()
        self.subgridvectorized = SubgridVectorized()
        self.subgriddb = SubgridDepthBased()
        if GPU:
            print('GPU (cupy) is not supported. numpy is used instead.')
            self.xp = numpy
        else:
            self.xp = numpy
        self.numHBasedTableLevels = numHBasedTableLevels

        self.readSubgridControlFile(subgridControlFilename)
        self.readMeshFile(self.control.meshFilename)

    ######################## Function to read subgrid control file ###################
    
    def readSubgridControlFile(self, subgridControlFilename):
        import re

        with open(subgridControlFilename) as ctrF:
            
            # skip a line
            ctrF.readline()

            # get output file name
            line = ctrF.readline().rstrip()
            line = re.split(' *= *',line)
            outputFilename = line[1]

            # get mesh filename
            line = ctrF.readline().rstrip()
            line = re.split(' *= *',line)
            meshFilename = line[1]

            # get list of elevation datasets
            line = ctrF.readline().rstrip()
            line = re.split(' *= *',line)
            numDEMs = int(line[1])
            demFilenameList = []
            for i in range(numDEMs):
                line = ctrF.readline().rstrip()
                demFilenameList.append(line)

            # get list of landcover datasets
            line = ctrF.readline().rstrip()
            line = re.split(' *= *',line)
            numLCs = int(line[1])
            landcoverFilenameList = []
            for i in range(numLCs):
                line = ctrF.readline().rstrip()
                landcoverFilenameList.append(line)

        self.control.outputFilename = outputFilename
        self.control.meshFilename = meshFilename
        self.control.numDEMs = numDEMs
        self.control.demFilenameList = demFilenameList
        self.control.numLCs = numLCs
        self.control.landcoverFilenameList = landcoverFilenameList
        self.control.loaded = True

    ##################### READ IN A MESH IN FORT.14 FORMAT #######################
    
    def readMeshFile(self, meshFilename):
        
        import pandas as pd
        import matplotlib.tri as mtri
        import numpy as np
        
        # initialize variables 
        
        x = []
        y = []
        z = []
        nodNum = []
        eleNum = []
        
        with open(meshFilename) as gridFile:
    
            gridFile.readline()
            line = gridFile.readline().split()
            numEle = int(line[0])
            numNod = int(line[1])
            
            triangles = []

            # import coordinates of points and elevations
            
            for i in range(numNod):
                
                line = gridFile.readline().split()
                nodNum.append(int(line[0]))
        
                x.append(float(line[1]))
                y.append(float(line[2]))
                z.append(float(line[3]))
        
            # import vertex order for each element
            # NOTE: the -1 in the triangles assembly is to make it 0 indexed
            
            for i in range(numEle):
                line = gridFile.readline().split()
                eleNum.append(int(line[0]))
                triangles.append([int(line[2])-1,int(line[3])-1,int(line[4])-1])
                                
        # data conversion
        triang = mtri.Triangulation(x,y,triangles)
                
        # put xyz into dataframe for ease of use
        gridXYZ = pd.DataFrame({'Node Number' : nodNum,
                     'Latitude' : y,
                     'Longitude' : x,
                     'Elevation' : z})

        self.mesh.coord = gridXYZ
        self.mesh.tri = triang.triangles
        self.mesh.numNod = numNod
        self.mesh.numEle = numEle
        self.mesh.loaded = True
        
    ############ CALCULATE AREA OF A TRIANGLE #######################################
    
    def triarea(x1,y1,x2,y2,x3,y3):
        
        area = abs((x1 * (y2-y3) + x2 * (y3 - y1)
                + x3 * (y1 -y2)) / 2.0)
    
        return area
    
    ############## CHECK IF A DEM CELL IS INSIDE A TRIANGULAR AREA ##################
    
    def isInside(x1,y1,x2,y2,x3,y3,x,y,difCriteria):
        from subgrid_calculator_dgp0 import SubgridCalculatorDGP0 as scm

        A = scm.triarea(x1,y1,x2,y2,x3,y3)
        A1 = scm.triarea(x,y,x2,y2,x3,y3)
        A2 = scm.triarea(x1,y1,x,y,x3,y3)
        A3 = scm.triarea(x1,y1,x2,y2,x,y)
        ADiff = abs(A - (A1 + A2 + A3))
    
        # this is a vectorized calulation that looks at all of the elements in 
        # an array
        mask = ADiff < difCriteria
    
        return mask
 
    #################### FUNCTION TO PROJECT WGS TO MERCATOR FOR CALCULATIONS #######################
    
    def projectMeshToMercator(lat, lon):
        
        import numpy as np
        
        r_major = 6378137.000
        x = r_major * np.radians(lon)
        scale = x/lon
        y = 180.0/np.pi * np.log(np.tan(np.pi/4.0 + lat * (np.pi/180.0)/2.0)) * scale
        return x, y

    ##################### FUNCTION TO IMPORT DEM #################################

    def importDEM(fileName):
        
        from osgeo import gdal
        import numpy as np
        
        # get dem info using gdal
        
        gdal_data = gdal.Open(fileName)
        gdal_band = gdal_data.GetRasterBand(1)
        nodataval = gdal_band.GetNoDataValue()
        demInfo = gdal_data.GetGeoTransform()
        xRes = demInfo[1] 
        x_upper_left = demInfo[0]
        y_upper_left = demInfo[3]
        yRes = demInfo[5]
        yDim = gdal_data.RasterYSize
        xDim = gdal_data.RasterXSize

        # convert to a numpy array
        
        z_array = gdal_data.ReadAsArray().astype(np.float64)
        
        # now create x,y coords array
        xCoords = np.zeros((xDim))
        yCoords = np.zeros((yDim))
        
        # gets points at the center of each raster cell
        
        for i in range(xDim):
            
            xCoords[i] = x_upper_left + (xRes)*i + (xRes/2)
        
        for i in range(yDim):
            
            yCoords[i] = y_upper_left + (yRes)*i + (yRes/2)
        
        # create a meshgrid (used for plotting)
        
        [X,Y] = np.meshgrid(xCoords,yCoords)
        
        # replace missing values if necessary
        # first set limits of zv
        
        if np.any(z_array == nodataval):
            z_array[z_array == nodataval] = np.nan
            
        return X,Y,z_array,xRes,yRes,xCoords,yCoords
    
    def importDEMv2(fileName):
        
        from osgeo import gdal
        import numpy as np
        
        # get dem info using gdal
        
        gdal_data = gdal.Open(fileName)
        gdal_band = gdal_data.GetRasterBand(1)
        nodataval = gdal_band.GetNoDataValue()
        demInfo = gdal_data.GetGeoTransform()
        xRes = demInfo[1] 
        x_upper_left = demInfo[0]
        y_upper_left = demInfo[3]
        yRes = demInfo[5]
        yDim = gdal_data.RasterYSize
        xDim = gdal_data.RasterXSize

        # convert to a numpy array
        
        z_array = gdal_data.ReadAsArray().astype(np.float32)
        
        # now create x,y coords array
        xCoords = np.zeros((xDim))
        yCoords = np.zeros((yDim))
        
        # gets points at the center of each raster cell
        
        for i in range(xDim):
            
            xCoords[i] = x_upper_left + (xRes)*i + (xRes/2)
        
        for i in range(yDim):
            
            yCoords[i] = y_upper_left + (yRes)*i + (yRes/2)
        
        # create a meshgrid (used for plotting)
        
        # [X,Y] = np.meshgrid(xCoords,yCoords)
        
        # replace missing values if necessary
        # first set limits of zv
        
        if np.any(z_array == nodataval):
            z_array[z_array == nodataval] = np.nan
            
        return z_array,xRes,yRes,xCoords,yCoords
    
    #################### FUNCTION TO PLOT SUBGRID VERTEX VARIABLES #####################

    def plotVertexSubgridVariable(meshObject,subgridVariable,levels=20):
        
        import matplotlib.pyplot as plt
        import cmocean
            
        fig1, ax1 = plt.subplots(figsize=(9,9))
        ax1.set_aspect('equal')
        tcf = ax1.tricontourf(meshObject[1], subgridVariable,cmap=cmocean.cm.rain,
                              levels=levels,extend='both')
        ax1.triplot(meshObject[1], color = 'k',linestyle='-',linewidth=0.25)
        cbar = fig1.colorbar(tcf,extendrect=True)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel('Elevation (m)', rotation=270,fontsize = 14)
        # ax1.set_title(title,fontsize=24)
        ax1.set_xlabel('Longitude',fontsize=20)
        ax1.set_ylabel('Latitude',fontsize=20)   
                    
    ########## CALCULATE SUBGRID CORRECTION FOR VECTORIZED STORAGE ##########

    def calculateLookupTableForVectorizedStorage(self):
        import sys
        import math
        import numpy as np
        import numpy as cp
        # import matplotlib.pyplot as plt
        # import matplotlib.style as style
        import time
        import matplotlib.pyplot as plt
        from subgrid_calculator_dgp0 import SubgridCalculatorDGP0 as scm

        # set the directory for everything to happen
        
        # inputDir = os.getcwd()
        
        startTot = time.perf_counter()
        
        # ignore any true divide errors
        np.seterr(invalid='ignore')
        
        # read in mesh
        numNod = self.mesh.numNod
        numEle = self.mesh.numEle

        # put elemental node coordinates in arrays
        xS = np.vstack((np.asarray(self.mesh.coord['Longitude'])[self.mesh.tri[:,0]],
                        np.asarray(self.mesh.coord['Longitude'])[self.mesh.tri[:,1]],
                        np.asarray(self.mesh.coord['Longitude'])[self.mesh.tri[:,2]])).T
        yS = np.vstack((np.asarray(self.mesh.coord['Latitude'])[self.mesh.tri[:,0]],
                        np.asarray(self.mesh.coord['Latitude'])[self.mesh.tri[:,1]],
                        np.asarray(self.mesh.coord['Latitude'])[self.mesh.tri[:,2]])).T
       
        # create a dictionary to hold the dem data
        elevationDict = {}
        
        for i in range(len(self.control.demFilenameList)):
            
            # read in DEM
            elevationData = scm.importDEM(self.control.demFilenameList[i])
            # x coordinates of DEM
            xDEMTemp = elevationData[0]
            # y coordinates of DEM
            yDEMTemp = elevationData[1]
            # elevations of DEM
            zDEMTemp = elevationData[2]
            # resolution in the x direction
            xDEMResTemp = elevationData[3]
            # resolution in the y direction
            yDEMResTemp = -1*elevationData[4]
            # x coordinates of DEM in 1D array
            xDEMCoordsTemp = elevationData[5]
            # y coordinates of DEM in 1D array
            yDEMCoordsTemp = elevationData[6]
            # deallocate DEM data 
            elevationData = None
            
            # get dem bounds to determine what dem to use when looping
            
            elevationDict["bounds%s"%i] = [np.min(xDEMCoordsTemp),
                                            np.max(xDEMCoordsTemp),
                                            np.min(yDEMCoordsTemp),
                                            np.max(yDEMCoordsTemp)]
            
            print('Finished reading DEM {0}.'.format(i))
            
        # deallocate arrays
        xDEMTemp = None
        yDEMTemp = None
        zDEMTemp = None
        xDEMResTemp = None
        yDEMResTemp = None
        xDEMCoordsTemp = None
        yDEMCoordsTemp = None
            
        # now find which elements are in which DEM's for processing later
        
        elementDict = {}
        
        totalEleInfoTable = np.empty((numEle,6))   # create empty array to hold dem values
        totalEleInfoTable[:,0] = np.arange(numEle) # polulate element numbers
        totalEleInfoTable[:,1] = np.min(xS,axis=1)  # populate minimum x values of the vertices of the element
        totalEleInfoTable[:,2] = np.max(xS,axis=1)  # populate maximum x values of the vertices of the element
        totalEleInfoTable[:,3] = np.min(yS,axis=1)  # populate minimum y values of the vertices of the element
        totalEleInfoTable[:,4] = np.max(yS,axis=1)  # populate maximum y values of the vertices of the element
        
        # loop through DEMs and determine which elements are in them
        # make sure to have DEMs in priority order meaning if you want to use fine
        # resolution in some areas have those dems listed first
        
        containedElementList0Index = []
        noElementDEMs = []
        
        for i in range(len(self.control.demFilenameList)):

            # find which elements are in the DEM
            totalEleInfoTable[:,5] = ((totalEleInfoTable[:,1]>elevationDict["bounds%s"%i][0])
                                      & ((totalEleInfoTable[:,2])<elevationDict["bounds%s"%i][1])
                                      & ((totalEleInfoTable[:,3])>elevationDict["bounds%s"%i][2])
                                      & ((totalEleInfoTable[:,4])<elevationDict["bounds%s"%i][3]))
        
            whichAreInside = list(np.where(totalEleInfoTable[:,5] == 1)[0])
            elementDict["DEM%s"%i] = totalEleInfoTable[whichAreInside,0].astype(int)        # store element numbers of the elements inside the DEM bound
            
            whichAreInsideActualEleNumber = totalEleInfoTable[whichAreInside,0].astype(int) # get the actual element numbers 
            totalEleInfoTable = np.delete(totalEleInfoTable,whichAreInside,axis=0)          # delete elements so we can skip those in the next dem

            # keep track if a dem does not have any elements inside to throw and exception later
            if len(whichAreInside) == 0:
                noElementDEMs.append(self.control.demFilenameList[i])
            
            containedElementList0Index.append(whichAreInsideActualEleNumber)    # create a list of elements within subgrid area
        
        # throw exception if a dem has no element and print those dem names
        if(len(noElementDEMs) != 0):
            for demName in noElementDEMs:
                print(demName)
            sys.exit('No elements in the above DEMs, throw those puppies out\n and their matching landcover!\n')
                        
        # now concatenate the lists from above
        containedElementList0Index = np.hstack(containedElementList0Index)
        
        # now delete double counted vertices and elements
        containedElementList0Index = np.unique(containedElementList0Index).astype(int)
        
        # now I want to create a list of 1s and 0s to show whether or not a
        # vertex or element is in the subgrid region
        
        binaryElementList = np.zeros(numEle)
        binaryElementList[containedElementList0Index] = 1
        
        # make an int array
        binaryElementList = binaryElementList.astype(int)
                
        # pre allocate arrays for subgrid quantities
        surfElevs = np.array([],dtype=np.float32)
        wetFraction = np.array([],dtype=np.float32)
        area = np.zeros((numEle)).astype(np.float32)
        totWatDepth = np.array([],dtype=np.float32)
        wetTotWatDepth = np.array([],dtype=np.float32)
        cf = np.array([],dtype=np.float32)
        minElevationEle = np.zeros(numEle).astype(np.float32)           # find lowest elevation in each element for use in variable phi
        maxElevationEle = np.zeros(numEle).astype(np.float32)           # find highest elevation in each element for use in variable phi

        # now fill the rows and columns of non subgrid vertices and elements with -99999
        area[np.where(binaryElementList == 0)] = -99999
        
        # fill min/max elevation
        minElevationEle[np.where(binaryElementList == 0)] = 99999
        maxElevationEle[np.where(binaryElementList == 0)] = -99999

        # allocate the index that maps element to location in vector
        elemIndex = np.zeros((numEle)).astype(np.int)
        elemNumLevel = np.zeros((numEle)).astype(np.int)
        elemIndexCnt = 0

        # these variables are used if you want level 1 corrections
        if self.level0andLevel1:
            rv = np.array([],dtype=np.float32)
            cmf = np.array([],dtype=np.float32)

        # specify buffer to use in area calculator
        # areaDif = 0.00001 
        areaDif = 0.00000001 
        
        # create variable to keep track of what DEM you have read in
        
        # first create a loop for DEMs
        for i in range(len(self.control.demFilenameList)):
            
            # reading in DEM again
            # all variables the same as before
            elevationData = scm.importDEM(self.control.demFilenameList[i])
            landcoverData = scm.importDEM(self.control.landcoverFilenameList[i])
            xDEM = elevationData[0]
            yDEM = elevationData[1]
            zDEM = elevationData[2]
            xDEMRes = elevationData[3]
            yDEMRes = -1*elevationData[4]
            xDEMCoords = elevationData[5]
            yDEMCoords = elevationData[6]
            elevationData = None # deallocate 
            nArray = landcoverData[2] # array of mannings n values
            landcoverData = None # deallocate
            # dictionary to translate between C-CAP and Manning's values
            # landCoverToManning = {0:0.02,2:0.12,3:0.12,4:0.12,5:0.035,6:0.1,7:0.05,8:0.035,
            #                   9:0.16,10:0.18,11:0.17,12:0.08,13:0.15,14:0.075,
            #                   15:0.06,16:0.15,17:0.07,18:0.05,19:0.03,20:0.03,
            #                   21:0.025,22:0.035,23:0.03,25:0.012}
            # change mannings conversion to match OM2D
            landCoverToManning = {0:0.02, 2:0.15, 3:0.10, 4:0.05, 5:0.02,
                                  6:0.037, 7:0.033, 8:0.034, 9:0.1, 10:0.11,
                                  11:0.1, 12:0.05, 13:0.1, 14:0.048, 15:0.045,
                                  16:0.1, 17:0.048, 18:0.045, 19:0.04,
                                  20:0.09, 21:0.02, 22:0.015, 23:0.015, 
                                  24:0.09, 25:0.01}
            # updated to SACS C-CAP mannings table
            # landCoverToManning = {0:0.025, 2:0.12, 3:0.10, 4:0.07, 5:0.035,
            #                       6:0.01, 7:0.055, 8:0.035, 9:0.16, 10:0.18,
            #                       11:0.17, 12:0.08, 13:0.15, 14:0.075, 15:0.07,
            #                       16:0.15, 17:0.07, 18:0.05, 19:0.03,
            #                       20:0.03, 21:0.025, 22:0.035, 23:0.03, 
            #                       24:0.09, 25:0.01}
            
            # landCoverValues = [0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25]
            # add landcover 24
            landCoverValues = [0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25]

            # convert values
            for value in landCoverValues:
                nArray[nArray == value] = landCoverToManning[value]
            
            
            # nArray[np.isnan(nArray)] = 0.012  # set nan values to 0.012 for ccap
            
            nArray[np.isnan(nArray)] = 0.02     # set nan values to 0.02
            
            # get a list of elements within this DEM
            elementList = elementDict["DEM%s"%i].astype(int)
            
            countElementLoop = 0

            # loop through the elements
            for ele in elementList:
            # for ele in elementList[0:100]:
                startTime = time.perf_counter()
            
                # get x and y points for vertices of this element
                currXPerimeterPoints = xS[ele,:]
                currYPerimeterPoints = yS[ele,:]
                
                countInElement = 0
                
                # get the bounds of the triangle
                xyArray = np.array((currXPerimeterPoints,currYPerimeterPoints)).T

                # convert to meters for area calculations
                xyCurrElementMeters = scm.projectMeshToMercator(xyArray[:,1],
                                                                xyArray[:,0])
                
                # calculate the area of the subsubelement
                eleArea = scm.triarea(xyCurrElementMeters[0][0],
                                      xyCurrElementMeters[1][0],
                                      xyCurrElementMeters[0][1],
                                      xyCurrElementMeters[1][1],
                                      xyCurrElementMeters[0][2],
                                      xyCurrElementMeters[1][2])
    
                # store element area to array
                area[ele] = eleArea
            
                # cut down DEM further to save in computation expense
                
                # max and min X and Y for the subsub element
                maxX = (np.max(currXPerimeterPoints) + 3*xDEMRes) # add a 3 cell buffer
                minX = (np.min(currXPerimeterPoints) - 3*xDEMRes)
                maxY = (np.max(currYPerimeterPoints) + 3*yDEMRes)
                minY = (np.min(currYPerimeterPoints) - 3*yDEMRes)
            
                # finds rows and column to cut down arrays
                minRow = np.where(yDEMCoords <= maxY)[0][0]
                maxRow = np.where(yDEMCoords >= minY)[0][-1]
                minCol = np.where(xDEMCoords >= minX)[0][0]
                maxCol = np.where(xDEMCoords <= maxX)[0][-1]
    
                # cut geotiff matrix down for faster processing
                xCutGeoTiffMatrix2 = xDEM[minRow:maxRow+1,minCol:maxCol+1]
                yCutGeoTiffMatrix2 = yDEM[minRow:maxRow+1,minCol:maxCol+1]
                zCutGeoTiffMatrix2 = zDEM[minRow:maxRow+1,minCol:maxCol+1]
                nCutGeoTiffMatrix2 = nArray[minRow:maxRow+1,minCol:maxCol+1]
            
                # for use in determining what cells lie within the element
                difCriteria = (np.ones((len(zCutGeoTiffMatrix2[:,0]),
                                len(zCutGeoTiffMatrix2[0,:])))*areaDif)
            
            
                # mask to find which cells are within an element
                mask = scm.isInside(currXPerimeterPoints[0],currYPerimeterPoints[0],
                                    currXPerimeterPoints[1],currYPerimeterPoints[1],
                                    currXPerimeterPoints[2],currYPerimeterPoints[2],
                                    xCutGeoTiffMatrix2,yCutGeoTiffMatrix2,
                                    difCriteria)
                
                # convert mask to cupy array
                mask = cp.asarray(mask)
                zCutGeoTiffMatrix2 = cp.asarray(zCutGeoTiffMatrix2)
                nCutGeoTiffMatrix2 = cp.asarray(nCutGeoTiffMatrix2)
                
                zCutGeoTiffMatrix2masked = zCutGeoTiffMatrix2[mask]
                nCutGeoTiffMatrix2masked = nCutGeoTiffMatrix2[mask]
                
                # get the min/max elevation of the element
                minElev = cp.min(zCutGeoTiffMatrix2masked)
                maxElev = cp.max(zCutGeoTiffMatrix2masked)

                # count how many cells are within subsubelement
                countIn = cp.count_nonzero(mask)
                
                # if there are no cells within the element the DEM is too coarse
                # you must decrease the DEM resolution in this area
                if countIn == 0:
                    sys.exit('DEM {0} resolution too coarse!'.format(i))
            
                # keep track of this for use later
                countInElement += countIn
            
                # create array of surface elevations
                # this is used to calculate the subgrid variables for varying water depths
                surfElevIncrement = 0.1
                minElevFloor = math.floor(minElev*10)/10
                surfaceElevations = np.arange(minElevFloor,maxElev+surfElevIncrement,surfElevIncrement)
                if len(surfaceElevations) < 2:
                    surfaceElevations = np.array([minElev,maxElev]).astype(np.float32)
                else:
                    surfaceElevations[0] = minElev
                    surfaceElevations[-1] = maxElev
                num_SfcElevs = len(surfaceElevations)   # number of surface elevations we are running
                surfaceElevations = cp.asarray(surfaceElevations)

                # preallocate arrays
                wetDryElementList = cp.zeros(num_SfcElevs)   # for wet area fraction phi
                totWatElementList = cp.zeros(num_SfcElevs)   # for grid total water depth H_G
                cfElementList = cp.zeros(num_SfcElevs)       # for coefficient of friction cf
                
                if self.level0andLevel1:
                    # for Level 1 calculation
                    rvBottomTermList = cp.zeros(num_SfcElevs)
    
                ################ BEGIN VECTORIZED CALCULATIONS ##############################
                # create a 3d surface array for use in calculations
                tempSurfaceElevArray = cp.ones((len(zCutGeoTiffMatrix2masked),
                                                    num_SfcElevs))*surfaceElevations

                # create a 3d manning array for use in calculations
            
                tempManningArray = nCutGeoTiffMatrix2masked[:,cp.newaxis]
            
                # subtract the bathymetry (2D Array) array from the surface 
                # elevations (3D array) to get total water depths over the 
                # element
                tempTotWatDepthArray = tempSurfaceElevArray -  zCutGeoTiffMatrix2masked[:,cp.newaxis]
                
                # find which of these cells are wet
                # add some tiny minimum water depth so we dont have cells with
                # like 10^-6 depths we will use 1 mm to start
                tempWetDryList = tempTotWatDepthArray > 0.001
                
                # count how many cells are wet
                tempCountWet = cp.count_nonzero(tempWetDryList,axis=0)

                ################ CALCULATE WET FRACTION #######################
                # keep track of how many cells are wet for use in averaging 
                # and wet fraction calulation
                wetDryElementList += tempCountWet

                ###############################################################

                #################### CALCULATING TOTAL WATER DEPTH ############
            
                # 0 out any dry cells
                tempTotWatDepthWetArray = tempTotWatDepthArray * tempWetDryList
                
                # integrate the wet total water depths for use in averaging later
                totWatElementList += cp.sum(tempTotWatDepthWetArray,axis=0)

                #################### CALCULATING MANNINGS n ###################
                # find the mannings for only wet areas then nan the rest for 
                # use in calculations 
                tempManningWetArray = tempManningArray * tempWetDryList

                ########### CALCULATE GRID AVERAGED CF FOR MANNING ############
                # calulate now for the element then sum later for use 
                # in other calulations
                tempcf = 9.81*tempManningWetArray**2/tempTotWatDepthWetArray**(1/3)

                ###############################################################
                
                if self.level0andLevel1:
                    ############ NOW CALCULATE RV FROM KENNEDY ET AL 2019 #########
                    # integrate only now and then calculate full rv later
                    # rvBottomTermList += np.nansum(tempTotWatDepthWetArray**(3/2)*\
                    #     (tempcf)**(-1/2),axis = 0)
                    rvBottomTermList += cp.nansum(tempTotWatDepthWetArray**(3/2)*\
                        (tempcf)**(-1/2),axis = 0)                    

                # finally sum and add cf for use in grid averaged calculation
                # tempcf = np.nansum(tempcf,axis=0)
                tempcf = cp.nansum(tempcf,axis=0)
                cfElementList += tempcf
            
                # Okay now we can finalize the values
                wetAvgTotWatDepth = totWatElementList/wetDryElementList
                gridAvgTotWatDepth = totWatElementList/countInElement
                wetFractionTemp = wetDryElementList/countInElement
                cfTemp = cfElementList/countInElement
                cmfTemp = (wetAvgTotWatDepth)*(wetAvgTotWatDepth/(rvBottomTermList/wetDryElementList))**2

                ####### convert to np array #######
                def get_ndarray(array):
                    if isinstance(array, np.ndarray):
                        return array
                    else:
                        return cp.ndarray.get(array)
                surfaceElevations = get_ndarray(surfaceElevations)
                gridAvgTotWatDepth = get_ndarray(gridAvgTotWatDepth)
                wetAvgTotWatDepth = get_ndarray(wetAvgTotWatDepth)
                wetFractionTemp = get_ndarray(wetFractionTemp)
                cfTemp = get_ndarray(cfTemp)
                if self.level0andLevel1:
                    cmfTemp = get_ndarray(cmfTemp)

                ####### store the values #######
                # store the index
                elemIndex[ele] = elemIndexCnt
                elemNumLevel[ele] = len(surfaceElevations)
                elemIndexCnt = elemIndexCnt + len(surfaceElevations)

                # store values
                surfElevs = np.append(surfElevs,surfaceElevations)
                totWatDepth = np.append(totWatDepth,gridAvgTotWatDepth)
                wetTotWatDepth = np.append(wetTotWatDepth,wetAvgTotWatDepth)
                wetFraction = np.append(wetFraction,wetFractionTemp)
                cf = np.append(cf,cfTemp)
                if self.level0andLevel1:
                    cmf = np.append(cmf,cmfTemp)

                minElevationEle[ele] = minElev
                maxElevationEle[ele] = maxElev
                
                countElementLoop += 1
                if countElementLoop%1000==0:
                    stopTime = time.perf_counter()
                    print("Finished Element {0} of {1} in DEM {2} took {3}".format(countElementLoop,len(elementList),i,stopTime - startTime))

        # add bottom limit on cf and cmf
        cf[cf<0.0025] = 0.0025
        cf[np.isnan(cf)] = 0.0025
        
        if self.level0andLevel1:
            cmf[cmf<0.0025] = 0.0025   
            cmf[np.isnan(cmf)] = 0.0025
            
        wetFraction[np.isnan(wetFraction)] = 0.0
        wetTotWatDepth[np.isnan(wetTotWatDepth)] = 0.0
        totWatDepth[np.isnan(totWatDepth)] = 0.0

        # sort the vectors in the ascending order of the elemental no.
        print('Sorting...')
        def sortarray(vals):
            sorted = np.zeros(len(vals)).astype(np.float32)
            sortedIndexStart = np.zeros(numEle).astype(np.int)
            sortedIndexEnd = np.zeros(numEle).astype(np.int)
            sortedIndexCnt = 0
            for ele in range(numEle):
                numLevels = elemNumLevel[ele]
                s = elemIndex[ele]
                e = elemIndex[ele]+numLevels
                snew = sortedIndexCnt
                enew = sortedIndexCnt+numLevels
                sorted[snew:enew] = vals[s:e]
                sortedIndexStart[ele] = snew
                sortedIndexEnd[ele] = enew
                sortedIndexCnt = enew
            return sorted, sortedIndexStart
        surfElevs, elemIndexSorted = sortarray(surfElevs)
        wetFraction, elemIndexSorted = sortarray(wetFraction)
        totWatDepth, elemIndexSorted = sortarray(totWatDepth)
        cfTemp, elemIndexSorted = sortarray(cfTemp)
        if self.level0andLevel1:
            cmfTemp, elemIndexSorted = sortarray(cmfTemp)
        print('Done.')

        # Find global min/max surface elevations
        minElevationGlobal = np.min(minElevationEle)
        maxElevationGlobal = np.max(maxElevationEle)

        # store the resulting values
        self.subgridvectorized.elemIndex = elemIndexSorted
        self.subgridvectorized.surfaceElevations = surfElevs
        self.subgridvectorized.wetFraction = wetFraction
        self.subgridvectorized.area = area
        self.subgridvectorized.totWatDepth = totWatDepth
        self.subgridvectorized.binaryElementList = binaryElementList
        self.subgridvectorized.minElevationEle = minElevationEle
        self.subgridvectorized.maxElevationEle = maxElevationEle
        self.subgridvectorized.minElevationGlobal = minElevationGlobal
        self.subgridvectorized.maxElevationGlobal = maxElevationGlobal
        self.subgridvectorized.cfElement = cf
        if self.level0andLevel1:
            self.subgridvectorized.cmfElement = cmf
        self.subgridvectorized.loaded = True
        
    ############## WRITE SUBGRID CORRECTION DATA TO NETCDF FOR VECTORIZED STORAGE ##############

    def writeSubgridLookupTableNetCDFForVectorizedStorage(self):
        from netCDF4 import Dataset
        import numpy as np

        # open netcdf file
        ncFile = Dataset(self.control.outputFilename, mode = 'w', format = 'NETCDF4')
        
        # create dimensions
        ncFile.createDimension('numEle',self.mesh.numEle)                               # element dimension
        ncFile.createDimension('numNode',self.mesh.numNod)                              # number of nodes in mesh
        ncFile.createDimension('vecLen',len(self.subgridvectorized.surfaceElevations))  # vector length
        ncFile.createDimension('global',1)                                              # globally single value variable
        
        # create variables
        # write variable dimensions transposed because FORTRAN will read them that way
        # only need to do it for the 3D arrays, handle 2 and 1D a different way.
        
        # indexes showing the locations of elements in the vectorized storage
        elemIndex = ncFile.createVariable('indexElementList',np.int,'numEle')
        
        # elemental wet area fraction
        wetFractionVarElement = ncFile.createVariable('wetFractionElement',np.float32,'vecLen')

        # elemental areas
        areaVar = ncFile.createVariable('area',np.float32,'numEle')
        
        # elemental grid averaged total water depth
        totWatDepthVar = ncFile.createVariable('totWatDepth',np.float32,'vecLen')
        
        # surface elevation array
        surfaceElevationsVar = ncFile.createVariable('surfaceElevations',np.float32,'vecLen')
        
        # elemental coefficient of friction level 0
        cfVarElement = ncFile.createVariable('cfElement',np.float32,'vecLen')

        # elemental coefficient of friction level 1
        if self.level0andLevel1:
            cmfVarElement = ncFile.createVariable('cmfElement',np.float32,'vecLen')

        # variables showing which elements are contained within the subgrid area
        binaryElementListVariable = ncFile.createVariable('binaryElementList',np.int,'numEle')
        
        # min/max Elevation
        minElevationEleVariable = ncFile.createVariable('minElevationElement',np.float32,'numEle')
        maxElevationEleVariable = ncFile.createVariable('maxElevationElement',np.float32,'numEle')
        minElevationGlobalVariable = ncFile.createVariable('minElevationGlobal',np.float32,'global')
        maxElevationGlobalVariable = ncFile.createVariable('maxElevationGlobal',np.float32,'global')
        
        elemIndex[:] = self.subgridvectorized.elemIndex
        wetFractionVarElement[:] = self.subgridvectorized.wetFraction
        areaVar[:] = self.subgridvectorized.area
        totWatDepthVar[:] = self.subgridvectorized.totWatDepth
        surfaceElevationsVar[:] = self.subgridvectorized.surfaceElevations
        cfVarElement[:] = self.subgridvectorized.cfElement
        if self.level0andLevel1:
            cmfVarElement[:] = self.subgridvectorized.cmfElement
        binaryElementListVariable[:] = self.subgridvectorized.binaryElementList
        minElevationEleVariable[:] = self.subgridvectorized.minElevationEle
        maxElevationEleVariable[:] = self.subgridvectorized.maxElevationEle
        minElevationGlobalVariable[:] = self.subgridvectorized.minElevationGlobal
        maxElevationGlobalVariable[:] = self.subgridvectorized.maxElevationGlobal
        
        ncFile.close()


    ######## READ SUBGRID CORRECTION DATA TO NETCDF FOR VECTORIZED STORAGE ########

    def readSubgridLookupTableNetCDFForVectorizedStorage(self):
        import numpy as np
        from netCDF4 import Dataset

        # open netcdf file
        lookupTable = Dataset(self.control.outputFilename)

        # read variables
        self.subgridvectorized.elemIndex = np.asarray(lookupTable['indexElementList'][:])
        self.subgridvectorized.surfaceElevations = np.asarray(lookupTable['surfaceElevations'][:])
        self.subgridvectorized.wetFraction = np.asarray(lookupTable['wetFractionElement'][:])
        self.subgridvectorized.area = np.asarray(lookupTable['area'][:])
        self.subgridvectorized.totWatDepth = np.asarray(lookupTable['totWatDepth'][:])
        self.subgridvectorized.cf = np.asarray(lookupTable['cfElement'][:]).T
        if self.level0andLevel1:
            self.subgridvectorized.cmf = np.asarray(lookupTable['cmfElement'][:]).T
        self.subgridvectorized.binaryElementList = np.asarray(lookupTable['binaryElementList'][:])
        self.subgridvectorized.minElevationEle = np.asarray(lookupTable['minElevationElement'][:])
        self.subgridvectorized.maxElevationEle = np.asarray(lookupTable['maxElevationElement'][:])
        self.subgridvectorized.minElevationGlobal = lookupTable['minElevationGlobal'][:]
        self.subgridvectorized.maxElevationGlobal = lookupTable['maxElevationGlobal'][:]

        lookupTable.close()

    ######## Plot subgrid correction data for vectorized storage ########

    def plot_subgrid_for_vectorizedstorage(self,imagePath):
        import numpy as np
        from matplotlib.patches import Circle, Wedge, Polygon
        from matplotlib.collections import PatchCollection
        import matplotlib.pyplot as plt

        # create polygons
        patches = []
        for i in range(self.mesh.numEle):
            t = self.mesh.tri[i,:]
            x = np.asarray(self.mesh.coord['Longitude'])[t]
            y = np.asarray(self.mesh.coord['Latitude'])[t]
            p = np.concatenate((x[:,np.newaxis],y[:,np.newaxis]), axis=1)

            polygon = Polygon(p, True)
            patches.append(polygon)

        def plotSingle(self,patches,ve,title,figfile):
            import numpy as np
            from matplotlib.patches import Circle, Wedge, Polygon
            from matplotlib.collections import PatchCollection
            import matplotlib.pyplot as plt

            p = PatchCollection(patches)
            p.set_array(ve)

            fig, ax = plt.subplots()
            ax.add_collection(p)
            fig.colorbar(p, ax=ax)
            plt.xlim([min(self.mesh.coord['Longitude']), max(self.mesh.coord['Longitude'])])
            plt.ylim([min(self.mesh.coord['Latitude']), max(self.mesh.coord['Latitude'])])
            plt.title(title)
            plt.show()
            fig.savefig("{:s}/sg_{:s}.png".format(imagePath,figfile),dpi=600,facecolor='white',edgecolor='none')

        def plotMulti(self,patches,numEle,elemIndex,surfaceElevs,minElevationEle,maxElevationEle,ve,surfElevForPlot,title,figfile):
            import numpy as np
            from matplotlib.patches import Circle, Wedge, Polygon
            from matplotlib.collections import PatchCollection
            import matplotlib.pyplot as plt

            def interpAt(elemIndex,surfaceElevs,minElevationEle,maxElevationEle,ve,se):
                import sys
                import math
                vinterp = np.zeros(numEle)
                for ele in range(numEle):
                    if se < minElevationEle[ele] or se > maxElevationEle[ele]:
                        vinterp[ele] = np.nan
                    else:
                        ista = elemIndex[ele]

                        if ele == numEle - 1:
                            iend = len(ve)
                        else:
                            iend = elemIndex[ele+1]
                        found = False
                        while True:
                            ii = math.floor((ista+iend)/2)
                            if surfaceElevs[ii] <= se and surfaceElevs[ii+1] >= se:
                                found = True
                                break
                            elif surfaceElevs[ii] > se:
                                iend = ii
                            elif surfaceElevs[ii+1] < se:
                                ista = ii
                        if found == False:
                            sys.exit("Matching slot not found")
                        r = (se - surfaceElevs[ii])/(surfaceElevs[ii+1] - surfaceElevs[ii])
                        vinterp[ele] = (1-r)*ve[ii] + r*ve[ii+1]
                return vinterp

            for i in range(0,len(surfElevForPlot)):
                se = surfElevForPlot[i] 
                vinterp = interpAt(elemIndex,surfaceElevs,minElevationEle,maxElevationEle,ve,se)

                p = PatchCollection(patches)
                p.set_array(vinterp)

                fig, ax = plt.subplots()
                ax.add_collection(p)
                fig.colorbar(p, ax=ax)
                plt.xlim([min(self.mesh.coord['Longitude']), max(self.mesh.coord['Longitude'])])
                plt.ylim([min(self.mesh.coord['Latitude']), max(self.mesh.coord['Latitude'])])
                plt.title("{:s} at Elevation {:.1f}".format(title,se))
                plt.show()
                fig.savefig("{:s}/sg_{:s}_{:03d}.png".format(imagePath,figfile,i),dpi=600,facecolor='white',edgecolor='none')

        surfElevForPlot = np.arange(-5.0,6.0,1.0)

        # mesh z
        vn = np.asarray(self.mesh.coord['Elevation'])
        ve = np.mean(vn[self.mesh.tri],axis=1)
        plotSingle(self,patches,ve,"Elemental Elevation","elevation")

        # area
        plotSingle(self,patches,self.subgridvectorized.area,"Elemental Area","area")

        # wetFraction
        plotMulti(self,patches,self.mesh.numEle,self.subgridvectorized.elemIndex,
                  self.subgridvectorized.surfaceElevations,
                  self.subgridvectorized.minElevationEle, self.subgridvectorized.maxElevationEle,
                  self.subgridvectorized.wetFraction,surfElevForPlot,"Elemental Wet Fraction","wetfraction")

        # totWatDepth
        plotMulti(self,patches,self.mesh.numEle,self.subgridvectorized.elemIndex,
                  self.subgridvectorized.surfaceElevations,
                  self.subgridvectorized.minElevationEle, self.subgridvectorized.maxElevationEle,
                  self.subgridvectorized.totWatDepth,surfElevForPlot,"Elemental Total Water Depth","totwatdepth")

        # binaryElementList
        plotSingle(self,patches,self.subgridvectorized.binaryElementList,"Binary Element List","binaryelemlist")

        # minElevationEle
        plotSingle(self,patches,self.subgridvectorized.minElevationEle,"Elemental Min Elevation","minelev")

        # maxElevationEle
        plotSingle(self,patches,self.subgridvectorized.maxElevationEle,"Elemental Max Elevation","maxelev")

        # cf
        plotMulti(self,patches,self.mesh.numEle,self.subgridvectorized.elemIndex,
                  self.subgridvectorized.surfaceElevations,
                  self.subgridvectorized.minElevationEle, self.subgridvectorized.maxElevationEle,
                  self.subgridvectorized.cf,surfElevForPlot,"Elemental Cf","cf")

        # cmf
        if self.level0andLevel1:
            plotMulti(self,patches,self.mesh.numEle,self.subgridvectorized.elemIndex,
                    self.subgridvectorized.surfaceElevations,
                    self.subgridvectorized.minElevationEle, self.subgridvectorized.maxElevationEle,
                    self.subgridvectorized.cmf,surfElevForPlot,"Elemental CMf","cmf")


