# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 16:17:03 2020

@author: jlwoodr3
"""

class subgridCalculatorDGP0main():

######################## Function to read subgrid control file ###################
    
    def readSubgridControlFile(subgridControlFilename):
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

        return outputFilename, meshFilename, numDEMs, demFilenameList, numLCs, landcoverFilenameList               

##################### READ IN A MESH IN FORT.14 FORMAT #######################
    
    def readMesh(self, meshFilename):
        
        import pandas as pd
        import matplotlib.tri as mtri
        import numpy as np
        
        # initialize variables 
        
        x = []
        y = []
        z = []
        vertNum = []
        eleNum = []
        triangles = []
        
        with open(meshFilename) as gridFile:
    
            gridFile.readline()
            line = gridFile.readline().split()
            numEle = int(line[0])
            numVert = int(line[1])
            
            # import coordinates of points and elevations
            
            for i in range(numVert):
                
                line = gridFile.readline().split()
                vertNum.append(int(line[0]))
        
                x.append(float(line[1]))
                y.append(float(line[2]))
                z.append(float(line[3]))
        
            # import vertex order for each element
            # NOTE: the -1 in the triangles assembly is to make it 0 indexed
            
            for i in range(numEle):
                
                line = gridFile.readline().split()
                eleNum.append(int(line[0]))
                triangles.append([int(line[2])-1,int(line[3])-1,int(line[4])-1])
                
            # triangulate mesh
            
        triang = mtri.Triangulation(x,y,triangles)
                
        # put xyz into dataframe for ease of use
        
        gridXYZ = pd.DataFrame({'Vertex Number' : vertNum,
                     'Latitude' : y,
                     'Longitude' : x,
                     'Elevation' : z})

        # self.gridXYZ = gridXYZ
        # self.gridElems = triang
        # self.numNode = numVert
        # self.numElem = numEle
            
        return gridXYZ, triang, numVert, numEle

############ CALCULATE AREA OF A TRIANGLE #######################################
    
    def triarea(x1,y1,x2,y2,x3,y3):
        
        area = abs((x1 * (y2-y3) + x2 * (y3 - y1)
                + x3 * (y1 -y2)) / 2.0)
    
        return area
    
############## CHECK IF A DEM CELL IS INSIDE A TRIANGULAR AREA ##################
    
    def isInside(x1,y1,x2,y2,x3,y3,x,y,difCriteria):
        from subgrid_calculator_dgp0 import subgridCalculatorDGP0main as scm

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
                    
########################## ADD SUBGRID CORRECTION WITH GPU ###################

    def calculateSubgridCorrectionwGPU(controlFilename):      

        import sys
        import numpy as np
        # import matplotlib.pyplot as plt
        # import matplotlib.style as style
        import time
        from netCDF4 import Dataset
        import cupy as cp
        import matplotlib.pyplot as plt
        from subgrid_calculator_dgp0 import subgridCalculatorDGP0main as scm

        # set the directory for everything to happen
        
        # inputDir = os.getcwd()
        
        startTot = time.perf_counter()
        
        # state if you want level0 and level1 corrections or just level 0
        
        level0andLevel1 = True # both if true, level0 only if false
        
        # add in control file
        
        # controlFile = 'E:/Fall2021/varible_dryPhi_tester/test_mesh_1/inputFiles/subgrid_control.txt'
        controlFile = controlFilename
        
        # ignore any true divide errors
        np.seterr(invalid='ignore')
        
        # read the control file
        (outputFilename, meshFilename, numDEMs,
        demFilenameList, numLCs, landcoverFilenameList) \
            = scm.readSubgridControlFile(controlFile)
                
        # read in mesh
        mesh = scm.readMesh(meshFilename)
        meshConnectivity = mesh[1].triangles
        numVert = mesh[2]
        numEle = mesh[3]

        # put elemental node coordinates in arrays
        xS = np.vstack((np.asarray(mesh[0]['Longitude'])[meshConnectivity[:,0]],
                        np.asarray(mesh[0]['Longitude'])[meshConnectivity[:,1]],
                        np.asarray(mesh[0]['Longitude'])[meshConnectivity[:,2]])).T
        yS = np.vstack((np.asarray(mesh[0]['Latitude'])[meshConnectivity[:,0]],
                        np.asarray(mesh[0]['Latitude'])[meshConnectivity[:,1]],
                        np.asarray(mesh[0]['Latitude'])[meshConnectivity[:,2]])).T
       
        # create a dictionary to hold the dem data
        elevationDict = {}
        
        for i in range(len(demFilenameList)):
            
            # read in DEM
            elevationData = scm.importDEM(demFilenameList[i])
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
        
        for i in range(len(demFilenameList)):

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
                noElementDEMs.append(demFilenameList[i])
            
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
        
        # create array of surface elevations
        # this is used to calculate the subgrid variables for varying water depths
        surfaceElevations = np.arange(-20,20.10,0.1)
        surfaceElevations = np.round(surfaceElevations,2)
        
        # now we start the computing 
        num_SfcElevs = len(surfaceElevations)   # number of surface elevations we are running
        nodesPerEle = 3                         # nodes per element
        
        # pre allocate arrays for subgrid quantities
        wetFraction = np.zeros((numEle,num_SfcElevs)).astype(np.float32)
        area = np.zeros((numEle)).astype(np.float32)
        totWatDepth = np.zeros((numEle,num_SfcElevs)).astype(np.float32)
        wetTotWatDepth = np.zeros((numEle,num_SfcElevs)).astype(np.float32)
        cf = np.zeros((numEle,num_SfcElevs)).astype(np.float32)
        maxElevationEle = np.zeros(numEle).astype(np.float32)           # find highest elevation in each element for use in variable phi
        

        # now fill the rows and columns of non subgrid vertices and elements with -99999
        wetFraction[np.where(binaryElementList == 0),:] = -99999
        area[np.where(binaryElementList == 0)] = -99999
        totWatDepth[np.where(binaryElementList == 0),:] = -99999
        wetTotWatDepth[np.where(binaryElementList == 0),:] = -99999
        cf[np.where(binaryElementList == 0),:] = -99999
        
        # fill max elevation
        maxElevationEle[np.where(binaryElementList == 0)] = -99999
        
        # these variables are used if you want level 1 corrections
        if level0andLevel1:
            rv = np.zeros((numEle,num_SfcElevs)).astype(np.float32)
            cmf = np.zeros((numEle,num_SfcElevs)).astype(np.float32)
            
            # now fill the rows and columns of non subgrid vertices and elements with -99999
            
            rv[np.where(binaryElementList == 0),:] = -99999
            cmf[np.where(binaryElementList == 0),:] = -99999

        # specify buffer to use in area calculator
        # areaDif = 0.00001 
        areaDif = 0.00000001 
        
        # create variable to keep track of what DEM you have read in
        
        # first create a loop for DEMs
        for i in range(len(demFilenameList)):
            
            # reading in DEM again
            # all variables the same as before
            elevationData = scm.importDEM(demFilenameList[i])
            landcoverData = scm.importDEM(landcoverFilenameList[i])
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
                startTime = time.perf_counter()
            
                # get x and y points for vertices of this element
                currXPerimeterPoints = xS[ele,:]
                currYPerimeterPoints = yS[ele,:]
                
                # preallocate arrays
                wetDryElementList = cp.zeros(num_SfcElevs)   # for wet area fraction phi
                totWatElementList = cp.zeros(num_SfcElevs)   # for grid total water depth H_G
                cfElementList = cp.zeros(num_SfcElevs)       # for coefficient of friction cf
                
                if level0andLevel1:
                    # for Level 1 calculation
                    rvBottomTermList = cp.zeros(num_SfcElevs)
    
                countInElement = 0
                
                # create a list of max elevation to keep track of elevation within 
                # each subelement then at the end of all the loops take the max
                # elevation and put it in the maxElevation array
                maxElevationTemp = []
            
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
                surfaceElevations = cp.asarray(surfaceElevations)
                
                zCutGeoTiffMatrix2masked = zCutGeoTiffMatrix2[mask]
                nCutGeoTiffMatrix2masked = nCutGeoTiffMatrix2[mask]
                
                # count how many cells are within subsubelement
                countIn = cp.count_nonzero(mask)
                
                # if there are no cells within the element the DEM is too coarse
                # you must decrease the DEM resolution in this area
                if countIn == 0:
                    sys.exit('DEM {0} resolution too coarse!'.format(i))
            
                # keep track of this for use later
                countInElement += countIn
            
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
                
                # max elevation
                maxElevationTemp.append(cp.max(zCutGeoTiffMatrix2masked))

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
                
                if level0andLevel1:
                    ############ NOW CALCULATE RV FROM KENNEDY ET AL 2019 #########
                    # integrate only now and then calculate full rv later
                    # rvBottomTermList += np.nansum(tempTotWatDepthWetArray**(3/2)*\
                    #     (tempcf)**(-1/2),axis = 0)
                    rvBottomTermList += cp.nansum(tempTotWatDepthWetArray**(3/2)*\
                        (tempcf)**(-1/2),axis = 0)                    

                # bottomFricEndTime = time.perf_counter()
                # print('Bottom friction time {} s'.format(bottomFricEndTime-bottomFricStartTime))
                # finally sum and add cf for use in grid averaged calculation
                # tempcf = np.nansum(tempcf,axis=0)
                tempcf = cp.nansum(tempcf,axis=0)
                cfElementList += tempcf
            
                # plot subarrea and dem cells inside
                # fig,ax = plt.subplots(1,1)
                # # ax.set_aspect('equal')
                # ax.plot(xCurrSubElement,yCurrSubElement)
                # masknumpy = cp.ndarray.get(mask)
                # img = ax.scatter(xCutGeoTiffMatrix2,yCutGeoTiffMatrix2,c=masknumpy)
                # fig.colorbar(img)
                # subsubCalcEndTime = time.perf_counter()
                # print(subsubCalcEndTime-subsubCalcStartTime)

                # Okay now we can finalize the values
                wetAvgTotWatDepth = totWatElementList/wetDryElementList
                gridAvgTotWatDepth = totWatElementList/countInElement
                wetFractionTemp = wetDryElementList/countInElement
                cfTemp = cfElementList/countInElement
                cmfTemp = (wetAvgTotWatDepth)*(wetAvgTotWatDepth/(rvBottomTermList/wetDryElementList))**2

                # give grid averaged total water depth
                totWatDepth[ele,:] = cp.ndarray.get(gridAvgTotWatDepth)

                # give wet averaged total water depth
                wetTotWatDepth[ele,:] = cp.ndarray.get(wetAvgTotWatDepth)

                # get wet area fraction for the subelement 
                wetFraction[ele,:] = cp.ndarray.get(wetFractionTemp)

                # get grid averaged coefficient of friction for the subelement
                ########### THIS IS WHAT I NORMALLY USE #################
                # convert back to numpy
                cf[ele,:] = cp.ndarray.get(cfTemp)

                # get the max elevation of the element
                maxElevationEle[ele] = max(maxElevationTemp)

                if level0andLevel1:
                    # get corrected bottom friction for level 1 corrections
                    cmf[ele,:] = cp.ndarray.get(cmfTemp)
                
                countElementLoop += 1
                if countElementLoop%1000==0:
                    stopTime = time.perf_counter()
                    print("Finished Element {0} of {1} in DEM {2} took {3}".format(countElementLoop,len(elementList),i,stopTime - startTime))
                
        # add bottom limit on cf and cmf
        cf[cf<0.0025] = 0.0025
        cf[np.isnan(cf)] = 0.0025
        
        if level0andLevel1:
            cmf[cmf<0.0025] = 0.0025   
            cmf[np.isnan(cmf)] = 0.0025
            
        wetFraction[np.isnan(wetFraction)] = 0.0
        wetTotWatDepth[np.isnan(wetTotWatDepth)] = 0.0
        totWatDepth[np.isnan(totWatDepth)] = 0.0

        # write netcdf file
        
        ncFile = Dataset(outputFilename, mode = 'w', format = 'NETCDF4')
        
        # create dimensions
        
        ncFile.createDimension('numEle',numEle) # element dimension
        ncFile.createDimension('numSfcElevs',len(surfaceElevations)) # number of surface elevations
        ncFile.createDimension('numNode',numVert) # number of nodes in mesh
        
        # create variables
        # write variable dimensions transposed because FORTRAN will read them that way
        # only need to do it for the 3D arrays, handle 2 and 1D a different way.
        
        # elemental wet area fraction
        wetFractionVarElement = ncFile.createVariable('wetFractionElement',np.float32,
                                                ('numSfcElevs','numEle'))
        # elemental areas
        areaVar = ncFile.createVariable('area',np.float32,('numEle'))
        
        # elemental grid averaged total water depth
        totWatDepthVar = ncFile.createVariable('totWatDepth',np.float32,
                                                ('numSfcElevs','numEle'))
        
        # surface elevation array
        surfaceElevationsVar = ncFile.createVariable('surfaceElevations',np.float32,'numSfcElevs')
        
        # variables showing which elements are contained within
        # the subgrid area
        
        binaryElementListVariable = ncFile.createVariable('binaryElementList',np.int,
                                                          ('numEle'))
        
        # write max Elevation
        maxElevationEleVariable = ncFile.createVariable('maxElevationElement',np.float32,
                                                      ('numEle'))
        
        if level0andLevel1:
            # elemental coefficient of friction level 1
            cmfVarElement = ncFile.createVariable('cmf',np.float32,
                                                 ('numEle','numSfcElevs'))
        
        wetFractionVarElement[:,:] = wetFraction.T
        areaVar[:] = area
        totWatDepthVar[:,:] = totWatDepth.T
        surfaceElevationsVar[:] = cp.ndarray.get(surfaceElevations)
        binaryElementListVariable[:] = binaryElementList
        # add max elevation cal
        maxElevationEleVariable[:] = maxElevationEle
        # cfVarElement[:,:,:] = cf
        if level0andLevel1:
            cmfVarElement[:,:] = cmf
        
        ncFile.close()
        
        endTot = time.perf_counter()
        
        print('Total Time = {} minutes'.format((endTot-startTot)/60))
