# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 16:17:03 2020

@author: jlwoodr3
"""

class subgridCalculatormain():

##################### READ IN A MESH IN FORT.14 FORMAT #######################
    
    def readMesh(meshFilename):
        
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

            
        return gridXYZ,triang, numVert, numEle

############ CALCULATE AREA OF A TRIANGLE #######################################
    
    def triarea(x1,y1,x2,y2,x3,y3):
        
        area = abs((x1 * (y2-y3) + x2 * (y3 - y1)
                + x3 * (y1 -y2)) / 2.0)
    
        return area
    
############## CHECK IF A DEM CELL IS INSIDE A TRIANGULAR AREA ##################
    
    def isInside(x1,y1,x2,y2,x3,y3,x,y,difCriteria):
    
        A = subgridCalculatormain.triarea(x1,y1,x2,y2,x3,y3)
        
        A1 = subgridCalculatormain.triarea(x,y,x2,y2,x3,y3)
    
        A2 = subgridCalculatormain.triarea(x1,y1,x,y,x3,y3)
    
        A3 = subgridCalculatormain.triarea(x1,y1,x2,y2,x,y)
    
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
                    
#################### FUNCTION TO PERFORM SUBGRID CALCULATIONS #######################    
    
    def calculateSubgridCorrection(controlFilename):
    
    
        import sys
        import numpy as np
        # import matplotlib.pyplot as plt
        # import matplotlib.style as style
        import time
        from netCDF4 import Dataset

        # set the directory for everything to happen
        
        # inputDir = os.getcwd()
        
        startTot = time.perf_counter()
        
        # state if you want level0 and level1 corrections or just level 0
        
        level0andLevel1 = True # both if true, level0 only if false
        
        # add in control file
        
        # controlFile = 'E:/Fall2021/varible_dryPhi_tester/test_mesh_1/inputFiles/subgrid_control.txt'
        controlFile = controlFilename
        
        # read the control file
        
        with open(controlFile) as ctrF:
            
            ctrF.readline()
            line = ctrF.readline().split()
            # get output file name
            outputFilename = line[2]
            line = ctrF.readline().split()
            # get mesh filename
            meshFilename = line[2]
            line = ctrF.readline().split()
            numDEMs = int(line[2])
            # get list of elevation datasets
            demFilenameList = []
            for i in range(numDEMs):
                line = ctrF.readline().split()
                demFilenameList.append(line[0])
            line = ctrF.readline().split()
            numLCs = int(line[2])
            # get list of landcover datasets
            landcoverFilenameList = []
            for i in range(numLCs):
                line = ctrF.readline().split()
                landcoverFilenameList.append(line[0])
                
        # read in mesh
        
        mesh = subgridCalculatormain.readMesh(meshFilename)
        
        # get connectivity table from mesh
        
        meshConnectivity = mesh[1].triangles
        
        # get number of vertices from mesh
        numVert = mesh[2]
        # get number of elements from mesh
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
            elevationData = subgridCalculatormain.importDEM(demFilenameList[i])
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
        
        # create empty array to hold dem values
        totalEleInfoTable = np.empty((mesh[3],6))
        # polulate element numbers
        totalEleInfoTable[:,0] = np.arange(mesh[3])
        # populate minimum x values of the vertices of the element
        totalEleInfoTable[:,1] = np.min(xS,axis=1)
        # populate maximum x values of the vertices of the element
        totalEleInfoTable[:,2] = np.max(xS,axis=1)
        # populate minimum y values of the vertices of the element
        totalEleInfoTable[:,3] = np.min(yS,axis=1)
        # populate maximum y values of the vertices of the element
        totalEleInfoTable[:,4] = np.max(yS,axis=1)
        
        # create a vertex info list
        
        totalVertInfoTable = np.empty((mesh[2],4))
        
        # populate element numbers
        
        totalVertInfoTable[:,0] = np.arange(mesh[2])
        
        # populate x values of vertices
        
        totalVertInfoTable[:,1] = mesh[0]['Longitude']
        
        # populate y values of vertices
        
        totalVertInfoTable[:,2] = mesh[0]['Latitude']
        
        # loop through DEMs and determine which elements are in them
        # make sure to have DEMs in priority order meaning if you want to use fine
        # resolution in some areas have those dems listed first
        
        containedElementList0Index = []
        containedVertexList0Index = []
        noElementDEMs = []
        
        for i in range(len(demFilenameList)):
        
            # fine which elements are in the DEM
            totalEleInfoTable[:,5] = ((totalEleInfoTable[:,1]>elevationDict["bounds%s"%i][0])
                                      & ((totalEleInfoTable[:,2])<elevationDict["bounds%s"%i][1])
                                      & ((totalEleInfoTable[:,3])>elevationDict["bounds%s"%i][2])
                                      & ((totalEleInfoTable[:,4])<elevationDict["bounds%s"%i][3]))
        
            whichAreInside = list(np.where(totalEleInfoTable[:,5] == 1)[0])
            elementDict["DEM%s"%i] = totalEleInfoTable[whichAreInside,0]
            # delete those elements from the total list
            # totalEleInfoTable = np.delete(totalEleInfoTable,whichAreInside,axis=0)
            # keep track if a dem does not have any elements inside to throw and exception later
            if len(whichAreInside) == 0:
                
                noElementDEMs.append(demFilenameList[i])
            # create a list of elements within subgrid area
            
            containedElementList0Index.append(whichAreInside)
            
            # create list of vertices within subgrid area
            
            totalVertInfoTable[:,3] = ((totalVertInfoTable[:,1]>elevationDict["bounds%s"%i][0])
                                      & ((totalVertInfoTable[:,1])<elevationDict["bounds%s"%i][1])
                                      & ((totalVertInfoTable[:,2])>elevationDict["bounds%s"%i][2])
                                      & ((totalVertInfoTable[:,2])<elevationDict["bounds%s"%i][3]))
            
            whichAreInside = list(np.where(totalVertInfoTable[:,3]==1)[0])
            
            containedVertexList0Index.append(whichAreInside)
            
            # keep track of vertices outside subgrid area
            
            # totalVertInfoTable = np.delete(totalVertInfoTable,whichAreInside,axis=0)
            
        # throw exception if a dem has no element and print those dem names
        
        if(len(noElementDEMs) != 0):
            
            for demName in noElementDEMs:
                
                print(demName)
        
            sys.exit('No elements in the above DEMs, throw those puppies out\n and their matching landcover!\n')
            
        # now concatenate the lists from above
        
        containedElementList0Index = np.hstack(containedElementList0Index)
        containedVertexList0Index = np.hstack(containedVertexList0Index)
        
        # now delete double counted vertices and elements
        
        containedElementList0Index = np.unique(containedElementList0Index).astype(int)
        containedVertexList0Index = np.unique(containedVertexList0Index).astype(int)
        
        # now I want to create a list of 1s and 0s to show whether or not a
        # vertex or element is in the subgrid region
        
        binaryElementList = np.zeros(mesh[3])
        binaryVertexList = np.zeros(mesh[2])
        
        binaryElementList[containedElementList0Index] = 1
        binaryVertexList[containedVertexList0Index] = 1
        
        # make an int array
        binaryElementList = binaryElementList.astype(int)
        binaryVertexList = binaryVertexList.astype(int)
            
        
        # find if there are any elements that do not have a DEM, or any that lie on a
        # boundary if so stop code
        
        # if(len(totalEleInfoTable[:,0]) != 0):
            
        #     sys.exit('Hold up! There are elements not totally inside any DEMs')
        
        
        # create array of surface elevations
        # this is used to calculate the subgrid variables for varying water depths
        
        surfaceElevations = np.arange(-20,20.10,0.1)
        surfaceElevations = np.round(surfaceElevations,2)
        
        # now we start the computing 
        
        # number of surface elevations we are running
        num_SfcElevs = len(surfaceElevations)
        
        # nodes per element
        nodesPerEle = 3
        
        # pre allocate arrays for subgrid quantities
        
        wetFraction = np.zeros((numEle,nodesPerEle,num_SfcElevs)).astype(np.float32)
        
        area = np.zeros((numEle,nodesPerEle)).astype(np.float32)
        
        totWatDepth = np.zeros((numEle,nodesPerEle,num_SfcElevs)).astype(np.float32)
        
        wetTotWatDepth = np.zeros((numEle,nodesPerEle,num_SfcElevs)).astype(np.float32)
        
        cf = np.zeros((numEle,nodesPerEle,num_SfcElevs)).astype(np.float32)
        
        maxElevationEle = np.zeros(numEle).astype(np.float32) # find highest elevation in each element for use in variable phi
        
        maxElevationSubEle = np.zeros((numEle,3)).astype(np.float32) # find the highest elevation in each subElement
        
        # now fill the rows and columns of non subgrid vertices and elements with -99999
        
        # wetFraction[totalEleInfoTable[:,0].astype(int),:,:] = -99999
        # area[totalEleInfoTable[:,0].astype(int),:] = -99999
        # totWatDepth[totalEleInfoTable[:,0].astype(int),:,:] = -99999
        # wetTotWatDepth[totalEleInfoTable[:,0].astype(int),:,:] = -99999
        # cf[totalEleInfoTable[:,0].astype(int),:,:] = -99999
        
        wetFraction[np.where(binaryElementList == 0),:,:] = -99999
        area[np.where(binaryElementList == 0)] = -99999
        totWatDepth[np.where(binaryElementList == 0),:,:] = -99999
        wetTotWatDepth[np.where(binaryElementList == 0),:,:] = -99999
        cf[np.where(binaryElementList == 0),:,:] = -99999
        
        # fill max elevation
        maxElevationEle[np.where(binaryElementList == 0)] = -99999
        maxElevationSubEle[np.where(binaryElementList == 0)] = -99999
        
        # these variables are used if you want level 1 corrections
        
        if level0andLevel1:
            
            rv = np.zeros((numEle,nodesPerEle,num_SfcElevs)).astype(np.float32)
            
            cmf = np.zeros((numEle,nodesPerEle,num_SfcElevs)).astype(np.float32)
            
            cadv = np.zeros((numEle,nodesPerEle,num_SfcElevs)).astype(np.float32)
            
            # now fill the rows and columns of non subgrid vertices and elements with -99999
            
            rv[np.where(binaryElementList == 0),:,:] = -99999
            cmf[np.where(binaryElementList == 0),:,:] = -99999
            cadv[np.where(binaryElementList == 0),:,:] = -99999
            
        

        # specify buffer to use in area calculator
        
        areaDif = 0.00001 
        
        # create variable to keep track of what DEM you have read in
        
        # first create a loop for DEMs
        
        for i in range(len(demFilenameList)):
            
            # reading in DEM again
            # all variables the same as before
            elevationData = subgridCalculatormain.importDEM(demFilenameList[i])
            landcoverData = subgridCalculatormain.importDEM(landcoverFilenameList[i])
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
            
            # landCoverValues = [0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25]
            # add landcover 24
            landCoverValues = [0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25]
            # covert values
        
            for value in landCoverValues:
            
                nArray[nArray == value] = landCoverToManning[value]
            
            # set nan values to 0.012 for ccap
        
            # nArray[np.isnan(nArray)] = 0.012
            
            # set nan values to 0.02
            
            nArray[np.isnan(nArray)] = 0.02
            
            # get a list of elements within this DEM
            # elementList = elementDict["DEM%s"%((len(demFilenameList)-1)-i)]
            elementList = elementDict["DEM%s"%i]
            
            countElementLoop = 0
            # loop through the elements
            for ele in elementList:
            # for ele in elementList[:2]:
                
                # cast ele to integer
                ele = int(ele)
                
                startTime = time.perf_counter()
            
                # get x and y points for vertices of this element
            
                currXPerimeterPoints = xS[ele,:]
                currYPerimeterPoints = yS[ele,:]
            
                # get x and y centroid
            
                currEleXCentroid = np.mean(currXPerimeterPoints)
                currEleYCenteroid = np.mean(currYPerimeterPoints)
                
                # loop through the 3 vertices and get averaged variables values within
                # each of the 3 sub elements
                
                for j in range(3):
                
                    # need to find mid points between vertices
                
                    nodeX = currXPerimeterPoints[j]
                    nodeY = currYPerimeterPoints[j]
                
                    otherXs = np.delete(currXPerimeterPoints,j)
                    otherYs = np.delete(currYPerimeterPoints,j)
        
                    midX1 = np.mean((nodeX,otherXs[0]))
                    midY1 = np.mean((nodeY,otherYs[0]))
                
                    midX2 = np.mean((nodeX,otherXs[1]))
                    midY2 = np.mean((nodeY,otherYs[1]))
        
                    # put the mid and nodal points in list so we have the perimeter of the
                    # sub areas
                
                    # xSubAreaPerimeter = [nodeX,midX1,currEleXCentroid,midX2]
                    # ySubAreaPerimeter = [nodeY,midY1,currEleYCenteroid,midY2]
                
                    # save this polygon for later use in patch plotter
                
                    # xyArray = np.array((xSubAreaPerimeter,ySubAreaPerimeter)).T
                
                    # polygon = Polygon(xyArray,True)
                    # patches.append(polygon)
                    
                    # get mid points of the sub element
                    midXs = [midX1,midX2]
                    midYs = [midY1,midY2]
                    
                    # now split that sub element into 2 sub sub elements to calcuate
                    # using triangles
                    numSubElements = 2
                
                    # preallocate arrays for use in subsubarea and subarea calculations
                    
                    # for wet area fraction phi
                    wetDrySubElementList = np.zeros(num_SfcElevs)
                    # for grid total water depth H_G
                    totWatSubElementList = np.zeros(num_SfcElevs)
                    # for coefficient of friction cf
                    cfSubElementList = np.zeros(num_SfcElevs)
                    # for wet total water depth H_W
                    # wetTotWatDepthSumList = np.zeros(num_SfcElevs)
                
                    if level0andLevel1:
                        # for Level 1 calculation
                        rvBottomTermList = np.zeros(num_SfcElevs)
                        
                        cadvMiddleTermList = np.zeros(num_SfcElevs)
        
                    countInSubElement = 0
                    
                    # create a list of max elevation to keep track of elevation within 
                    # each subelement then at the end of all the loops take the max
                    # elevation and put it in the maxElevation array
                
                    maxElevationTemp = []
                
                    # loop through each subsubelement then add quantities together
                    # to get averages over entire subelement
                
                    for k in range(numSubElements):
                    
                        # get the bounds of the triangle for each subsubArea
                    
                        xCurrSubElement = [nodeX, currEleXCentroid, midXs[k]]
                        yCurrSubElement = [nodeY, currEleYCenteroid, midYs[k]]
                        
                        xySubArray = np.array((xCurrSubElement,yCurrSubElement)).T
                        
                        # convert to meters for area calculations
                        
                        xyCurrSubElementMeters = subgridCalculatormain.projectMeshToMercator(xySubArray[:,1],
                                                                          xySubArray[:,0])
                    
                        
                        # calculate the area of the subsubelement
                        
                        subsubEleArea = subgridCalculatormain.triarea(xyCurrSubElementMeters[0][0],
                                                                      xyCurrSubElementMeters[1][0],
                                                                      xyCurrSubElementMeters[0][1],
                                                                      xyCurrSubElementMeters[1][1],
                                                                      xyCurrSubElementMeters[0][2],
                                                                      xyCurrSubElementMeters[1][2])
                    
                    
                        # add area subelement area calculation to area array
                    
                        area[ele,j] += subsubEleArea
                    
                        # cut down DEM further to save in computation expense
                        
                        # max and min X and Y for the subsub element
                        maxX = (np.max(xCurrSubElement) + 3*xDEMRes) # add a 3 cell buffer
                        minX = (np.min(xCurrSubElement) - 3*xDEMRes)
                    
                        maxY = (np.max(yCurrSubElement) + 3*yDEMRes)
                        minY = (np.min(yCurrSubElement) - 3*yDEMRes)
                    
        
                      # finds rows and column to cut down arrays
                    
                        minRow = np.where(yDEMCoords <= maxY)[0][0]
                        maxRow = np.where(yDEMCoords >= minY)[0][-1]
            
                        minCol = np.where(xDEMCoords >= minX)[0][0]
                        maxCol = np.where(xDEMCoords <= maxX)[0][-1]
                    
            
                        # cut geotiff matrix down for faster processing
                    
                        xCutGeoTiffMatrix2 = xDEM[minRow:maxRow+1,minCol:maxCol+1]#.astype(np.float32)
                        yCutGeoTiffMatrix2 = yDEM[minRow:maxRow+1,minCol:maxCol+1]#.astype(np.float32)
                        zCutGeoTiffMatrix2 = zDEM[minRow:maxRow+1,minCol:maxCol+1]#.astype(np.float32)
                        nCutGeoTiffMatrix2 = nArray[minRow:maxRow+1,minCol:maxCol+1]#.astype(np.float32)
                    
        
                        # for use in determining what cells lie within the subsubelement
                        difCriteria = (np.ones((len(zCutGeoTiffMatrix2[:,0]),
                                            len(zCutGeoTiffMatrix2[0,:])))*areaDif)
                    
                    
                        # mask to find which cells are within a subsubelement
                    
                        mask = subgridCalculatormain.isInside(xCurrSubElement[0],yCurrSubElement[0],
                                                              xCurrSubElement[1],yCurrSubElement[1],
                                                              xCurrSubElement[2],yCurrSubElement[2],
                                                              xCutGeoTiffMatrix2,yCutGeoTiffMatrix2,
                                                              difCriteria)
                    
                        # count how many cells are within subsubelement
                    
                        countIn = np.count_nonzero(mask)
                        
                        # if there are no cells within the element the DEM is too coarse
                        # you must decrease the DEM resolution in this area
                        if countIn == 0:
                        
                            sys.exit('DEM {0} resolution too coarse!'.format(i))
                    
                        # keep track of this for use later
                    
                        countInSubElement += countIn
                    
        ################ BEGIN VECTORIZED CALCULATIONS ##############################
                    
                        # create a 3d surface array for use in calculations
                    
                        tempSurfaceElevArray = (np.ones((len(zCutGeoTiffMatrix2[:,0]),
                                                      len(zCutGeoTiffMatrix2[0,:]),
                                                      num_SfcElevs))*surfaceElevations)
                    
                        # create a 3d manning array for use in calculations
                    
                        tempManningArray = nCutGeoTiffMatrix2[:,:,np.newaxis]
                    
                        # subtract the bathymetry (2D Array) array from the surface 
                        # elevations (3D array) to get total water depths over the 
                        # sub element 
                    
                        tempTotWatDepthArray = np.subtract(tempSurfaceElevArray,zCutGeoTiffMatrix2[:,:,None])
                        
                        # max elevation
                        
                        maxElevationTemp.append(np.max(zCutGeoTiffMatrix2[mask]))
                    
                        # only look at cells within the subsubelement
                    
                        tempTotWatDepthArray = tempTotWatDepthArray[mask]
                        tempManningArray = tempManningArray[mask]
        
                        # find which of these cells are wet
                    
                        # tempWetDryList = tempTotWatDepthArray > 0.0
                        # add some tiny minimum water depth so we dont have cells with
                        # like 10^-6 depths we will use 1 mm to start
                        tempWetDryList = tempTotWatDepthArray > 0.001
                    
                        # count how many cells are wet
                    
                        tempCountWet = np.count_nonzero(tempWetDryList,axis=0)
                    
                    
                    ################ CALCULATE WET FRACTION #######################
                    
                        # keep track of how many cells are wet for use in averaging 
                        # and wet fraction calulation
                    
                        wetDrySubElementList += tempCountWet
                    
                    ###############################################################
        
                    #################### CALCULATING TOTAL WATER DEPTH ############
                    
                        # nan out any dry cells
                        tempTotWatDepthWetArray = tempTotWatDepthArray * tempWetDryList
                        
                        # integrate the wet total water depths for use in averaging later
                        totWatSubElementList += np.nansum(tempTotWatDepthWetArray,axis=0)
                        
                        #################### CALCULATING MANNINGS n ###################
                        
                        # find the mannings for only wet areas then nan the rest for 
                        # use in calculations 
                        
                        tempManningWetArray = tempManningArray * tempWetDryList
                        
                        tempManningWetArray[tempManningWetArray == 0.0] = np.nan
                        
                        ########### CALCULATE GRID AVERAGED CF FOR MANNING ############
                        
                        # calulate now for the subsubelement then sum later for use 
                        # in other calulations
                        tempcf = 9.81*tempManningWetArray**2/tempTotWatDepthWetArray**(1/3)
        
                        ###############################################################
                        
                        if level0andLevel1:
                            ############ NOW CALCULATE RV FROM KENNEDY ET AL 2019 #########
            
                            # integrate only now and then calculate full rv later
                            
                            rvBottomTermList += np.nansum(tempTotWatDepthWetArray**(3/2)*\
                                (tempcf)**(-1/2),axis = 0)
            
                            ########### CALCULATE Part of C ADVECTION ####################
            
                            cadvMiddleTermList += np.nansum(tempTotWatDepthWetArray**2/tempcf
                                                            ,axis=0)
                    
                        
                        # finally sum and add cf for use in grid averaged calculation
                        tempcf = np.nansum(tempcf,axis=0)
                        cfSubElementList += tempcf
                    
                    
                    # Okay now we can calculate the averaged variable values for each
                    # of the 3 polygonal sub elements within a real adcirc element
                    
                    # go ahead and average wet and grid total water depth for use in calulations
                    wetDrySubElementList[wetDrySubElementList == 0] = np.nan
                    wetAvgTotWatDepth = totWatSubElementList/wetDrySubElementList
                    gridAvgTotWatDepth = totWatSubElementList/countInSubElement
                    
                    # give grid averaged total water depth
                    totWatDepth[ele,j,:] = gridAvgTotWatDepth
                    
                    # give wet averaged total water depth
                    wetTotWatDepth[ele,j,:] = wetAvgTotWatDepth
                    
                    # get wet area fraction for the subelement 
                    wetFraction[ele,j,:] = wetDrySubElementList/countInSubElement
                    
                    # get grid averaged coefficient of friction for the subelement
                    ######### JUST TESTING THIS #################
                    # cf[i,j,:] = cfSubElementList/countInSubElement
                    
                    # actually on inspection of implementation, I think I need to solve
                    # for wet averaged cf
                    ########### THIS IS WHAT I NORMALLY USE #################
                    cf[ele,j,:] = cfSubElementList/wetDrySubElementList
                
                    # get the max elevation of the subgelement
                    maxElevationSubEle[ele,j] = max(maxElevationTemp)
                    
                    if level0andLevel1:
                        # for this rv we average the bottom and top terms separatly over the whole subelement
                        rv[ele,j,:] = wetAvgTotWatDepth/(rvBottomTermList/wetDrySubElementList)
                        
                        # get corrected bottom friction for level 1 corrections
                        # cmf[i,j,:] = (gridAvgTotWatDepth)*rv[i,j,:]**2 # this is incorrect from the Kennedy et al 2019 paper
                        cmf[ele,j,:] = (wetAvgTotWatDepth)*rv[ele,j,:]**2
        
                        # get advection correction for level 1 corrections
                        cadv[ele,j,:] = (1/(wetAvgTotWatDepth)) \
                            *(cadvMiddleTermList/wetDrySubElementList)*rv[ele,j,:]**2
                        
                #### ADD MAX ELEVATION ###
                    
                maxElevationEle[ele] = np.max(maxElevationSubEle[ele,:])
                
                countElementLoop += 1
                stopTime = time.perf_counter()
                print("Finished Element {0} of {1} in DEM {2} took {3}".format(countElementLoop,len(elementList),i,stopTime - startTime))
            
                
        # add bottom limit on cf, cmf. and cadv
        cf[cf<0.0025] = 0.0025
        cf[np.isnan(cf)] = 0.0025
        
        if level0andLevel1:
            cmf[cmf<0.0025] = 0.0025   
            cmf[np.isnan(cmf)] = 0.0025
            cadv[np.isnan(cadv)] = 1.0    
            
        wetFraction[np.isnan(wetFraction)] = 0.0
        wetTotWatDepth[np.isnan(wetTotWatDepth)] = 0.0
        totWatDepth[np.isnan(totWatDepth)] = 0.0
        
        
        # put elemental quantities on vertices
        
        wetFractionVertex = np.zeros((numVert,num_SfcElevs))
        wetTotWatDepthVertex = np.zeros((numVert,num_SfcElevs))
        gridTotWatDepthVertex = np.zeros((numVert,num_SfcElevs))
        vertexArea = np.zeros((numVert,num_SfcElevs))
        cfVertex = np.zeros((numVert,num_SfcElevs))
        
        # add max elevation
        maxElevationVertex = np.ones(numVert)*-99999
    
        
        # now fill non subgrid spaces with -99999
        
        wetFractionVertex[np.where(binaryVertexList == 0),:] = -99999
        wetTotWatDepthVertex[np.where(binaryVertexList == 0),:] = -99999
        gridTotWatDepthVertex[np.where(binaryVertexList == 0),:] = -99999
        vertexArea[np.where(binaryVertexList == 0),:] = -99999
        cfVertex[np.where(binaryVertexList == 0),:] = -99999
        maxElevationVertex[np.where(binaryVertexList == 0)] = -99999

        
        if level0andLevel1:
            cmfVertex = np.zeros((numVert,num_SfcElevs))
            cadvVertex = np.zeros((numVert,num_SfcElevs))
            
            # now fill non subgrid spaces with -99999
            
            cmfVertex[np.where(binaryVertexList == 0),:] = -99999
            cadvVertex[np.where(binaryVertexList == 0),:] = -99999
        
        # loop through elements and sum there vertex quantities    
        for j in range(numEle):
            
            if(binaryElementList[j] == 1):
                
                nm1 = meshConnectivity[j,0]
                nm2 = meshConnectivity[j,1]
                nm3 = meshConnectivity[j,2]
                
                phi1 = wetFraction[j,0,:]
                phi2 = wetFraction[j,1,:]
                phi3 = wetFraction[j,2,:]
                
                HW1 = wetTotWatDepth[j,0,:]
                HW2 = wetTotWatDepth[j,1,:]
                HW3 = wetTotWatDepth[j,2,:]
                
                HG1 = totWatDepth[j,0,:]
                HG2 = totWatDepth[j,1,:]
                HG3 = totWatDepth[j,2,:]
                
                cf1 = cf[j,0,:]
                cf2 = cf[j,1,:]
                cf3 = cf[j,2,:]
                
                # add max elevation
                if(maxElevationVertex[nm1] < maxElevationSubEle[j,0]):
                    
                    maxElevationVertex[nm1] = maxElevationSubEle[j,0]
            
                if(maxElevationVertex[nm2] < maxElevationSubEle[j,1]):
                    
                    maxElevationVertex[nm2] = maxElevationSubEle[j,1]
                    
                if(maxElevationVertex[nm3] < maxElevationSubEle[j,2]):
                    
                    maxElevationVertex[nm3] = maxElevationSubEle[j,2]
                
                vertexArea[nm1,:] += area[j,0]
                vertexArea[nm2,:] += area[j,1]
                vertexArea[nm3,:] += area[j,2]
                
                # sum element quantities on vertex 
                # multiply by area to get area weighted average
                
                wetFractionVertex[nm1,:] += phi1 * area[j,0]
                wetFractionVertex[nm2,:] += phi2 * area[j,1]
                wetFractionVertex[nm3,:] += phi3 * area[j,2]
                
                wetTotWatDepthVertex[nm1,:] += HW1 * area[j,0]
                wetTotWatDepthVertex[nm2,:] += HW2 * area[j,1]
                wetTotWatDepthVertex[nm3,:] += HW3 * area[j,2]
                
                gridTotWatDepthVertex[nm1,:] += HG1 * area[j,0]
                gridTotWatDepthVertex[nm2,:] += HG2 * area[j,1]
                gridTotWatDepthVertex[nm3,:] += HG3 * area[j,2]
                
                cfVertex[nm1,:] += cf1 * area[j,0]
                cfVertex[nm2,:] += cf2 * area[j,1]
                cfVertex[nm3,:] += cf3 * area[j,2]
            
                if level0andLevel1:
                
                    cmf1 = cmf[j,0,:]
                    cmf2 = cmf[j,1,:]
                    cmf3 = cmf[j,2,:]
                    
                    cadv1 = cadv[j,0,:]
                    cadv2 = cadv[j,1,:]
                    cadv3 = cadv[j,2,:]
                
                    cmfVertex[nm1,:] += cmf1 * area[j,0]
                    cmfVertex[nm2,:] += cmf2 * area[j,1]
                    cmfVertex[nm3,:] += cmf3 * area[j,2]
                
                    cadvVertex[nm1,:] += cadv1 * area[j,0]
                    cadvVertex[nm2,:] += cadv2 * area[j,1]
                    cadvVertex[nm3,:] += cadv3 * area[j,2]
        
        # now average all of these by the vertex areas
        wetFractionVertex[np.where(binaryVertexList == 1)] = wetFractionVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        wetTotWatDepthVertex[np.where(binaryVertexList == 1)] = wetTotWatDepthVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        gridTotWatDepthVertex[np.where(binaryVertexList == 1)] = gridTotWatDepthVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        cfVertex[np.where(binaryVertexList == 1)] = cfVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        # cfwetVertex = cfwetVertex/vertexArea
        if level0andLevel1:
            
            cmfVertex[np.where(binaryVertexList == 1)] = cmfVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
            # cmfgridVertex = cmfgridVertex/vertexArea
            cadvVertex[np.where(binaryVertexList == 1)] = cadvVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
            # cadvgridVertex = cadvgridVertex/vertexArea
        
        # now we need to check to see if there are any vertices in the subgrid
        # but not connected to any elements in the subgrid 
        
        for i in range(numVert):
            
            # for vertices that are in the subgrid but may or may not be
            # connected to any subgrid elements
            if(binaryVertexList[i]==1):
                # find connected elements
                check = np.where(meshConnectivity == i)[0]
                # check to see if any of the elements are in the subgrid
                inSubgrid = any(binaryElementList[check] == 1)
                
                if(inSubgrid != True):
                    
                    # change the vertex to not be in the subgrid
                    binaryVertexList[i] = 0
                    # change the vertex variables to be -99999
                    wetFractionVertex[i,:] = -99999
                    wetTotWatDepthVertex[i,:] = -99999
                    gridTotWatDepthVertex[i,:] = -99999
                    cfVertex[i,:] = -99999
                    if(level0andLevel1):
                        cmfVertex[i,:] = -99999
                        cadvVertex[i,:] = -99999
                        
            # otherwise if vertex is outside subgrid area
            # just make everything -99999
            else:
                
                # change the vertex variables to be -99999
                wetFractionVertex[i,:] = -99999
                wetTotWatDepthVertex[i,:] = -99999
                gridTotWatDepthVertex[i,:] = -99999
                cfVertex[i,:] = -99999
                if(level0andLevel1):
                    cmfVertex[i,:] = -99999
                    cadvVertex[i,:] = -99999

        # write netcdf file
        
        # if level0andLevel1:
        #     outputFilename = 'E:/Fall2021/varible_dryPhi_tester/test_mesh_1/inputFiles/avgVars_01_increment.nc'
        # else:
        #     outputFilename = 'E:/Fall2021/varible_dryPhi_tester/test_mesh_1/inputFiles/avgVarslevel0only_01_increment.nc'
        
        ncFile = Dataset(outputFilename, mode = 'w', format = 'NETCDF4')
        
        # create dimensions
        
        ncFile.createDimension('numEle',numEle) # element dimension
        ncFile.createDimension('numVert',3) # number of vertices per element
        ncFile.createDimension('numSfcElevs',len(surfaceElevations)) # number of surface elevations
        ncFile.createDimension('numNode',numVert) # number of nodes in mesh
        
        # create variables
        # write variable dimensions transposed because FORTRAN will read them that way
        # only need to do it for the 3D arrays, handle 2 and 1D a different way.
        
        # vertex wet area fraction
        wetFractionVarVertex = ncFile.createVariable('wetFractionVertex',np.float32,
                                            ('numNode','numSfcElevs'))
        # elemental wet area fraction
        wetFractionVarElement = ncFile.createVariable('wetFractionElement',np.float32,
                                               ('numSfcElevs','numVert','numEle'))
        # elemental areas
        areaVar = ncFile.createVariable('area',np.float32,('numEle','numVert'))
        
        # elemental grid averaged total water depth
        totWatDepthVar = ncFile.createVariable('totWatDepth',np.float32,
                                               ('numSfcElevs','numVert','numEle'))
        
        # # vertex minimum wet water depths
        # phiMinDepthVar = ncFile.createVariable('phiMinDepth',np.float32,'numNode')
        
        # surface elevation array
        surfaceElevationsVar = ncFile.createVariable('surfaceElevations',np.float32,'numSfcElevs')
        
        # vertex wet total water depth
        wetTotWatDepthVarVertex = ncFile.createVariable('wetTotWatDepthVertex',np.float32,
                                                         ('numNode','numSfcElevs'))
        
        # vertex grid total water depth
        gridTotWatDepthVarVertex = ncFile.createVariable('gridTotWatDepthVertex',np.float32,
                                                         ('numNode','numSfcElevs'))
        
        # vertex coefficient of friction level 0
        cfVarVertex = ncFile.createVariable('cfVertex',np.float32,
                                            ('numNode','numSfcElevs'))
        
        # variables showing which elements and vertices are contained within
        # the subgrid area
        
        binaryElementListVariable = ncFile.createVariable('binaryElementList',np.int,
                                                          ('numEle'))
        binaryVertexListVariable = ncFile.createVariable('binaryVertexList',np.int,
                                                          ('numNode'))
        
        # write max Elevation
        
        maxElevationEleVariable = ncFile.createVariable('maxElevationElement',np.float32,
                                                     ('numEle'))
        maxElevationVertexVariable = ncFile.createVariable('maxElevationVertex',
                                                           np.float32,
                                                           ('numNode'))
        
        if level0andLevel1:
            # vertex coefficient of friction level 1
            cmfVarVertex = ncFile.createVariable('cmfVertex',np.float32,
                                                 ('numNode','numSfcElevs'))
            
            # elemental advection correction
            cadvVar = ncFile.createVariable('cadv',np.float32,
                                            ('numSfcElevs','numVert','numEle'))
        
        
        wetFractionVarVertex[:,:] = wetFractionVertex
        wetFractionVarElement[:,:,:] = wetFraction.T
        areaVar[:,:] = area
        totWatDepthVar[:,:,:] = totWatDepth.T
        # phiMinDepthVar[:] = phiMinDepth
        surfaceElevationsVar[:] = surfaceElevations
        wetTotWatDepthVarVertex[:,:] = wetTotWatDepthVertex
        gridTotWatDepthVarVertex[:,:] = gridTotWatDepthVertex
        cfVarVertex[:,:] = cfVertex
        binaryElementListVariable[:] = binaryElementList
        binaryVertexListVariable[:] = binaryVertexList
        # add max elevation cal
        maxElevationEleVariable[:] = maxElevationEle
        maxElevationVertexVariable[:] = maxElevationVertex
        # cfVarElement[:,:,:] = cf
        if level0andLevel1:
            
            cmfVarVertex[:,:] = cmfVertex
            # cmfVarElement[:,:,:] = cmf
            cadvVar[:,:,:] = cadv.T
        
        ncFile.close()
        
        endTot = time.perf_counter()
        
        print('Total Time = {} minutes'.format((endTot-startTot)/60))
        
        # return phi1011Check, phi1057Check, phi1056Check, phi1013Check, phi1010Check, phi1053Check
        
        # # try using a patch plot
         
        # plot vertex quantities
        
        # fig,ax = plt.subplots(3,1,figsize=(12,9))
        
        # ax[0].plot(surfaceElevations,cfVertex[186,:],label='<Cf>W')
        # # ax[0].plot(surfaceElevations,cfwetVertex[490,:],label='<Cf>W')
        # if level0andLevel1:
        #     ax[0].plot(surfaceElevations,cmfVertex[186,:],label='<H>W*Rv^2')
        #     # ax[0].plot(surfaceElevations,cmfgridVertex[490,:],label='<H>G*Rv^2')
        #     ax[0].set_title('cf vs cmf')
        # else:
        #     ax[0].set_title('cf')
            
        # ax[0].legend()
        
        # ax[1].plot(surfaceElevations,wetFractionVertex[186,:])     
        # ax[1].set_title('phi')   
        
        # if level0andLevel1:
        #     ax[2].plot(surfaceElevations,cadvVertex[186,:],label='wetavg')   
        #     # ax[2].plot(surfaceElevations,cadvgridVertex[490,:],label='gridavg middleG')   
        #     ax[2].set_title('cadv')  
        #     ax[2].set_ylim((0,2))
            
        # cmap = cmocean.cm.rain
        # colorLevels = np.arange(0,1.05,0.05)
        # cbarTicks = np.arange(0,1.1,0.1)
            
        # adcircToolsModule.resultsPlotter.timestep63plotter(wetFractionVertex[:,200],mesh[1],cmap,
        #                               colorLevels,cbarTicks,'phi Vertex')
        
################ FUNCTION TO DECREASE THE SIZE OF THE LOOKUP TABLE ############

    def simplifyLookupTable(lookupTableFilename, outputFilename):
        
        # this script will try and reduce the size of lookup tables 
        
        import numpy as np
        import netCDF4 as nc
        
        # read in the lookup table
        
        lookupTable = nc.Dataset(lookupTableFilename)
        
        # import variables
        
        wetAreaFractionElement = np.asarray(lookupTable['wetFractionElement'][:])
        totalWaterDepthEle = np.asarray(lookupTable['totWatDepth'][:])
        cadv = np.asarray(lookupTable['cadv'])
        depths = np.asarray(lookupTable['surfaceElevations'][:])
        wetAreaFractionVertex = np.asarray(lookupTable['wetFractionVertex'])
        gridTotWatDepthVertex = np.asarray(lookupTable['gridTotWatDepthVertex'])
        wetTotWatDepthVertex = np.asarray(lookupTable['wetTotWatDepthVertex'])
        area = np.asarray(lookupTable['area'])
        cfVertex = np.asarray(lookupTable['cfVertex'])
        cmfVertex = np.asarray(lookupTable['cmfVertex'])
        binaryElementList = np.asarray(lookupTable['binaryElementList'][:])
        binaryVertexList = np.asarray(lookupTable['binaryVertexList'][:])
        maxElevationEle = np.asarray(lookupTable['maxElevationElement'][:])
        maxElevationVertex = np.asarray(lookupTable['maxElevationVertex'][:])
        
        # create desired phi list
        
        desiredPhiList = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        
        # create empty arrays for the reduced tables
        
        depthsEleForLookup = np.zeros((11,3,len(wetAreaFractionElement[0,0,:])))
        HEleForLookup = np.zeros((11,3,len(wetAreaFractionElement[0,0,:])))
        cadvForLookup = np.zeros((11,3,len(wetAreaFractionElement[0,0,:])))

        # now loop through the elements 
        
        for i in range(len(wetAreaFractionElement[0,0,:])):
            
            for j in range(3):
                
                currPhiArray = wetAreaFractionElement[:,j,i]
                
                # check to see if the element is always fully wet
                
                if(all(currPhiArray == 1.0)):
                    
                    # set all to really low depths so that in the code you know to set at 1
                
                    depthsEleForLookup[:,j,i] = np.min(depths)
                    
                    # set the total water depth to the first value so that we have a reference
                    # point to add zeta to
                    
                    HEleForLookup[:,j,i] = totalWaterDepthEle[0,j,i]
                    
                    # set cadvection to the last value in the lookup table
                    
                    cadvForLookup[:,j,i] = cadv[-1,j,i]
                    
                # check to see if the element is always fully dry
                    
                elif(all(currPhiArray == 0.0)):
                    
                    # set all to really high dpeths so that in code you know to set to 0.0
                    
                    depthsEleForLookup[:,j,i] = np.max(depths) 
                    HEleForLookup[:,j,i] = 0.0
                    
                    # set cadv to the first value of the lookup table
                    
                    cadvForLookup[:,j,i] = cadv[0,j,i]
        
                # if the element is partially wet for all values
                            
                elif(all(currPhiArray > 0.0) and all(currPhiArray < 1.0)):
                    
                    # set to always dry
                    
                    depthsEleForLookup[:,j,i] = np.max(depths) 
                    HEleForLookup[:,j,i] = 0.0
                    cadvForLookup[:,j,i] = cadv[0,j,i]
                            
                # if partially wet to fully wet
                
                elif(all(currPhiArray > 0.0) and any(currPhiArray == 1.0)):
                    
                    # set to always wet
                    
                    depthsEleForLookup[:,j,i] = np.min(depths)
                    
                    # set the total water depth to the first value so that we have a reference
                    # point to add zeta to
                    
                    HEleForLookup[:,j,i] = totalWaterDepthEle[0,j,i]
                    
                    # set cadvection to the last value in the lookup table
                    
                    cadvForLookup[:,j,i] = cadv[-1,j,i]
                            
                # if fully dry to partially wet
                
                elif(all(currPhiArray < 1.0) and any(currPhiArray == 0.0)):
                    
                    # set to always dry
                    
                    depthsEleForLookup[:,j,i] = np.max(depths)  
                    HEleForLookup[:,j,i] = 0.0
                    cadvForLookup[:,j,i] = cadv[0,j,i]
                    
                # classical example where an element starts dry (phi == 0)
                # and becomes wet (phi == 1)
                    
                elif(any(currPhiArray == 0) and any(currPhiArray == 1)):
                
                    for k in range(len(desiredPhiList)):
                    
                        desiredPhi = desiredPhiList[k]
                        
                        # for phi == 0.0 you want exactly where that is in the currPhiArray
                        
                        if(desiredPhi == 0.0):
                            
                            equalTo = np.where(currPhiArray == desiredPhi)[0][-1]
                            depthsEleForLookup[k,j,i] = depths[equalTo]
                            HEleForLookup[k,j,i] = totalWaterDepthEle[equalTo,j,i]
                            cadvForLookup[k,j,i] = cadv[equalTo,j,i]
                            
                        # for phi == 1.0 you want exactly where that is in the currPhiArray
                            
                        elif(desiredPhi == 1.0):
                            
                            equalTo = np.where(currPhiArray == desiredPhi)[0][0]
                            depthsEleForLookup[k,j,i] = depths[equalTo]
                            HEleForLookup[k,j,i] = totalWaterDepthEle[equalTo,j,i]
                            cadvForLookup[k,j,i] = cadv[equalTo,j,i]
                            
                        # now look for the points in between 0 and 1
                            
                        else:
                            
                            greaterThan = np.where(currPhiArray > desiredPhi)[0][0]
                            lessThan = np.where(currPhiArray < desiredPhi)[0][-1]
                            
                            depthsEleForLookup[k,j,i] = (((desiredPhi - currPhiArray[lessThan])
                                                         /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                         *(depths[greaterThan] - depths[lessThan])
                                                         + (depths[lessThan]))
                            
                            HEleForLookup[k,j,i] = (((desiredPhi - currPhiArray[lessThan])
                                                         /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                         *(totalWaterDepthEle[greaterThan,j,i] - totalWaterDepthEle[lessThan,j,i])
                                                         + (totalWaterDepthEle[lessThan,j,i]))
                            
                            cadvForLookup[k,j,i] = (((desiredPhi - currPhiArray[lessThan])
                                                         /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                         *(cadv[greaterThan,j,i] - cadv[lessThan,j,i])
                                                         + (cadv[lessThan,j,i]))
                            
                            
                                     
            print('Element {0} completed'.format(i))
            
        # fill the elements that are not contained in the subgrid region with -99999
        depthsEleForLookup[:,:,np.where(binaryElementList==0)[0]] = -99999
        HEleForLookup[:,:,np.where(binaryElementList==0)[0]] = -99999
        cadvForLookup[:,:,np.where(binaryElementList==0)[0]] = -99999
        
        ncFile = nc.Dataset(outputFilename, mode = 'w', format = 'NETCDF4')
        
        # create dimensions
        
        numEle = len(cadv[0,0,:])
        numVert = len(gridTotWatDepthVertex[:,0])
        level0andLevel1 = True
        
        ncFile.createDimension('numEle',numEle) # element dimension
        ncFile.createDimension('numVert',3) # number of vertices per element
        ncFile.createDimension('numPhi',11) # number of possible phi values
        ncFile.createDimension('numSfcElevs',len(depths)) # number of surface elevations
        ncFile.createDimension('numNode',numVert) # number of nodes in mesh
        
        # create variables
        # write variable dimensions transposed because FORTRAN will read them that way
        # only need to do it for the 3D arrays, handle 2 and 1D a different way.
        
        # wetAreaFraction set
        
        phiSet = ncFile.createVariable('phiSet',np.float32,
                                       'numPhi')
        # vertex wet area fraction
        wetFractionVarVertex = ncFile.createVariable('wetFractionVertex',np.float32,
                                            ('numNode','numSfcElevs'))
        # elemental wet area fraction
        # wetFractionVarElement = ncFile.createVariable('wetFractionElement',np.float32,
        #                                        ('numPhi','numEle'))
        
        wetFractionVarDepths = ncFile.createVariable('wetFractionDepths',np.float32,
                                               ('numPhi','numVert','numEle'))
        # elemental areas
        areaVar = ncFile.createVariable('area',np.float32,('numEle','numVert'))
        
        # elemental grid averaged total water depth
        totWatDepthVar = ncFile.createVariable('totWatDepth',np.float32,
                                               ('numPhi','numVert','numEle'))
        
        # # vertex minimum wet water depths
        # phiMinDepthVar = ncFile.createVariable('phiMinDepth',np.float32,'numNode')
        
        # surface elevation array
        surfaceElevationsVar = ncFile.createVariable('surfaceElevations',np.float32,'numSfcElevs')
        
        # vertex wet total water depth
        wetTotWatDepthVarVertex = ncFile.createVariable('wetTotWatDepthVertex',np.float32,
                                                         ('numNode','numSfcElevs'))
        
        # vertex grid total water depth
        gridTotWatDepthVarVertex = ncFile.createVariable('gridTotWatDepthVertex',np.float32,
                                                         ('numNode','numSfcElevs'))
        
        # vertex coefficient of friction level 0
        cfVarVertex = ncFile.createVariable('cfVertex',np.float32,
                                            ('numNode','numSfcElevs'))
        
        # variables showing which elements and vertices are contained within
        # the subgrid area
                
        binaryElementListVariable = ncFile.createVariable('binaryElementList',np.int,
                                                                  ('numEle'))
        binaryVertexListVariable = ncFile.createVariable('binaryVertexList',np.int,
                                                                  ('numNode'))
        
        # write max Elevation
        
        maxElevationEleVariable = ncFile.createVariable('maxElevationElement',np.float32,
                                                     ('numEle'))
        
        maxElevationVertexVariable = ncFile.createVariable('maxElevationVertex',np.float32,
                                                     ('numNode'))
        
        if level0andLevel1:
            # vertex coefficient of friction level 1
            cmfVarVertex = ncFile.createVariable('cmfVertex',np.float32,
                                                 ('numNode','numSfcElevs'))
            
            # elemental advection correction
            cadvVar = ncFile.createVariable('cadv',np.float32,
                                            ('numPhi','numVert','numEle'))
        
        
        phiSet[:] = desiredPhiList
        wetFractionVarVertex[:,:] = wetAreaFractionVertex
        wetFractionVarDepths[:,:,:] = depthsEleForLookup
        areaVar[:,:] = area
        totWatDepthVar[:,:,:] = HEleForLookup
        # phiMinDepthVar[:] = phiMinDepth
        surfaceElevationsVar[:] = depths
        wetTotWatDepthVarVertex[:,:] = wetTotWatDepthVertex
        gridTotWatDepthVarVertex[:,:] = gridTotWatDepthVertex
        cfVarVertex[:,:] = cfVertex
        binaryElementListVariable[:] = binaryElementList
        binaryVertexListVariable[:] = binaryVertexList
        # add max elevation cal
        maxElevationEleVariable[:] = maxElevationEle
        maxElevationVertexVariable[:] = maxElevationVertex
        
        
        # cfVarElement[:,:,:] = cf
        if level0andLevel1:
            
            cmfVarVertex[:,:] = cmfVertex
            # cmfVarElement[:,:,:] = cmf
            cadvVar[:,:,:] = cadvForLookup
        
        ncFile.close()
                
########################## ADD SUBGRID CORRECTION WITH GPU ###################

    def calculateSubgridCorrectionwGPU(controlFilename):      

        import sys
        import numpy as np
        # import matplotlib.pyplot as plt
        # import matplotlib.style as style
        import time
        from netCDF4 import Dataset
        import cupy as cp

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
        
        with open(controlFile) as ctrF:
            
            ctrF.readline()
            line = ctrF.readline().split()
            # get output file name
            outputFilename = line[2]
            line = ctrF.readline().split()
            # get mesh filename
            meshFilename = line[2]
            line = ctrF.readline().split()
            numDEMs = int(line[2])
            # get list of elevation datasets
            demFilenameList = []
            for i in range(numDEMs):
                line = ctrF.readline().split()
                demFilenameList.append(line[0])
            line = ctrF.readline().split()
            numLCs = int(line[2])
            # get list of landcover datasets
            landcoverFilenameList = []
            for i in range(numLCs):
                line = ctrF.readline().split()
                landcoverFilenameList.append(line[0])
                
        # read in mesh
        
        mesh = subgridCalculatormain.readMesh(meshFilename)
        
        # get connectivity table from mesh
        
        meshConnectivity = mesh[1].triangles
        
        # get number of vertices from mesh
        numVert = mesh[2]
        # get number of elements from mesh
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
            elevationData = subgridCalculatormain.importDEM(demFilenameList[i])
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
        
        # create empty array to hold dem values
        totalEleInfoTable = np.empty((mesh[3],6))
        # polulate element numbers
        totalEleInfoTable[:,0] = np.arange(mesh[3])
        # populate minimum x values of the vertices of the element
        totalEleInfoTable[:,1] = np.min(xS,axis=1)
        # populate maximum x values of the vertices of the element
        totalEleInfoTable[:,2] = np.max(xS,axis=1)
        # populate minimum y values of the vertices of the element
        totalEleInfoTable[:,3] = np.min(yS,axis=1)
        # populate maximum y values of the vertices of the element
        totalEleInfoTable[:,4] = np.max(yS,axis=1)
        
        # create a vertex info list
        
        totalVertInfoTable = np.empty((mesh[2],4))
        
        # populate element numbers
        
        totalVertInfoTable[:,0] = np.arange(mesh[2])
        
        # populate x values of vertices
        
        totalVertInfoTable[:,1] = mesh[0]['Longitude']
        
        # populate y values of vertices
        
        totalVertInfoTable[:,2] = mesh[0]['Latitude']
        
        # loop through DEMs and determine which elements are in them
        # make sure to have DEMs in priority order meaning if you want to use fine
        # resolution in some areas have those dems listed first
        
        containedElementList0Index = []
        containedVertexList0Index = []
        noElementDEMs = []
        
        for i in range(len(demFilenameList)):
        
            # fine which elements are in the DEM
            totalEleInfoTable[:,5] = ((totalEleInfoTable[:,1]>elevationDict["bounds%s"%i][0])
                                      & ((totalEleInfoTable[:,2])<elevationDict["bounds%s"%i][1])
                                      & ((totalEleInfoTable[:,3])>elevationDict["bounds%s"%i][2])
                                      & ((totalEleInfoTable[:,4])<elevationDict["bounds%s"%i][3]))
        
            whichAreInside = list(np.where(totalEleInfoTable[:,5] == 1)[0])
            elementDict["DEM%s"%i] = totalEleInfoTable[whichAreInside,0]
            # delete those elements from the total list
            # totalEleInfoTable = np.delete(totalEleInfoTable,whichAreInside,axis=0)
            # keep track if a dem does not have any elements inside to throw and exception later
            if len(whichAreInside) == 0:
                
                noElementDEMs.append(demFilenameList[i])
            
            # create a list of elements within subgrid area
            
            containedElementList0Index.append(whichAreInside)
            
            # create list of vertices within subgrid area
            
            totalVertInfoTable[:,3] = ((totalVertInfoTable[:,1]>elevationDict["bounds%s"%i][0])
                                      & ((totalVertInfoTable[:,1])<elevationDict["bounds%s"%i][1])
                                      & ((totalVertInfoTable[:,2])>elevationDict["bounds%s"%i][2])
                                      & ((totalVertInfoTable[:,2])<elevationDict["bounds%s"%i][3]))
            
            whichAreInside = list(np.where(totalVertInfoTable[:,3]==1)[0])
            
            containedVertexList0Index.append(whichAreInside)
            
            # keep track of vertices outside subgrid area
            
            # totalVertInfoTable = np.delete(totalVertInfoTable,whichAreInside,axis=0)
        
        # throw exception if a dem has no element and print those dem names
        
        if(len(noElementDEMs) != 0):
            
            for demName in noElementDEMs:
                
                print(demName)
        
            sys.exit('No elements in the above DEMs, throw those puppies out\n and their matching landcover!\n')
                        
        # now concatenate the lists from above
        
        containedElementList0Index = np.hstack(containedElementList0Index)
        containedVertexList0Index = np.hstack(containedVertexList0Index)
        
        # now delete double counted vertices and elements
        
        containedElementList0Index = np.unique(containedElementList0Index).astype(int)
        containedVertexList0Index = np.unique(containedVertexList0Index).astype(int)
        
        # now I want to create a list of 1s and 0s to show whether or not a
        # vertex or element is in the subgrid region
        
        binaryElementList = np.zeros(mesh[3])
        binaryVertexList = np.zeros(mesh[2])
        
        binaryElementList[containedElementList0Index] = 1
        binaryVertexList[containedVertexList0Index] = 1
        
        # make an int array
        binaryElementList = binaryElementList.astype(int)
        binaryVertexList = binaryVertexList.astype(int)
            
        
        # find if there are any elements that do not have a DEM, or any that lie on a
        # boundary if so stop code
        
        # if(len(totalEleInfoTable[:,0]) != 0):
            
        #     sys.exit('Hold up! There are elements not totally inside any DEMs')
        
        
        # create array of surface elevations
        # this is used to calculate the subgrid variables for varying water depths
        
        surfaceElevations = np.arange(-20,20.10,0.1)
        surfaceElevations = np.round(surfaceElevations,2)
        
        # now we start the computing 
        
        # number of surface elevations we are running
        num_SfcElevs = len(surfaceElevations)
        
        # nodes per element
        nodesPerEle = 3
        
        # pre allocate arrays for subgrid quantities
        
        wetFraction = np.zeros((numEle,nodesPerEle,num_SfcElevs)).astype(np.float32)
        
        area = np.zeros((numEle,nodesPerEle)).astype(np.float32)
        
        totWatDepth = np.zeros((numEle,nodesPerEle,num_SfcElevs)).astype(np.float32)
        
        wetTotWatDepth = np.zeros((numEle,nodesPerEle,num_SfcElevs)).astype(np.float32)
        
        cf = np.zeros((numEle,nodesPerEle,num_SfcElevs)).astype(np.float32)
        
        maxElevationEle = np.zeros(numEle).astype(np.float32) # find highest elevation in each element for use in variable phi
        
        maxElevationSubEle = np.zeros((numEle,3)).astype(np.float32) # find the highest elevation in each subElement
        

        # now fill the rows and columns of non subgrid vertices and elements with -99999
        
        # wetFraction[totalEleInfoTable[:,0].astype(int),:,:] = -99999
        # area[totalEleInfoTable[:,0].astype(int),:] = -99999
        # totWatDepth[totalEleInfoTable[:,0].astype(int),:,:] = -99999
        # wetTotWatDepth[totalEleInfoTable[:,0].astype(int),:,:] = -99999
        # cf[totalEleInfoTable[:,0].astype(int),:,:] = -99999
        
        wetFraction[np.where(binaryElementList == 0),:,:] = -99999
        area[np.where(binaryElementList == 0)] = -99999
        totWatDepth[np.where(binaryElementList == 0),:,:] = -99999
        wetTotWatDepth[np.where(binaryElementList == 0),:,:] = -99999
        cf[np.where(binaryElementList == 0),:,:] = -99999
        
        # fill max elevation
        maxElevationEle[np.where(binaryElementList == 0)] = -99999
        maxElevationSubEle[np.where(binaryElementList == 0)] = -99999
        
        # these variables are used if you want level 1 corrections
        
        if level0andLevel1:
            
            rv = np.zeros((numEle,nodesPerEle,num_SfcElevs)).astype(np.float32)
            
            cmf = np.zeros((numEle,nodesPerEle,num_SfcElevs)).astype(np.float32)
            
            cadv = np.zeros((numEle,nodesPerEle,num_SfcElevs)).astype(np.float32)
            
            # now fill the rows and columns of non subgrid vertices and elements with -99999
            
            rv[np.where(binaryElementList == 0),:,:] = -99999
            cmf[np.where(binaryElementList == 0),:,:] = -99999
            cadv[np.where(binaryElementList == 0),:,:] = -99999
            
        

        # specify buffer to use in area calculator
        
        areaDif = 0.00001 
        
        # create variable to keep track of what DEM you have read in
        
        # first create a loop for DEMs
        
        for i in range(len(demFilenameList)):
            
            # reading in DEM again
            # all variables the same as before
            elevationData = subgridCalculatormain.importDEM(demFilenameList[i])
            landcoverData = subgridCalculatormain.importDEM(landcoverFilenameList[i])
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
            
            # landCoverValues = [0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25]
            # add landcover 24
            landCoverValues = [0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25]
            # covert values
        
            for value in landCoverValues:
            
                nArray[nArray == value] = landCoverToManning[value]
            
            # set nan values to 0.012 for ccap
        
            # nArray[np.isnan(nArray)] = 0.012
            
            # set nan values to 0.02
            
            nArray[np.isnan(nArray)] = 0.02
            
            # get a list of elements within this DEM
            # elementList = elementDict["DEM%s"%((len(demFilenameList)-1)-i)]
            elementList = elementDict["DEM%s"%i].astype(int)
            
            countElementLoop = 0
            # loop through the elements
            for ele in elementList:
            # for i in range(1):
            # for ele in elementList[:2]:
                
                # cast ele to integer
                # ele = int(elementList[i])
                
                startTime = time.perf_counter()
            
                # get x and y points for vertices of this element
            
                currXPerimeterPoints = xS[ele,:]
                currYPerimeterPoints = yS[ele,:]
            
                # get x and y centroid
            
                currEleXCentroid = np.mean(currXPerimeterPoints)
                currEleYCenteroid = np.mean(currYPerimeterPoints)
                
                # loop through the 3 vertices and get averaged variables values within
                # each of the 3 sub elements
                
                for j in range(3):
                
                    # need to find mid points between vertices
                
                    nodeX = currXPerimeterPoints[j]
                    nodeY = currYPerimeterPoints[j]
                
                    otherXs = np.delete(currXPerimeterPoints,j)
                    otherYs = np.delete(currYPerimeterPoints,j)
        
                    midX1 = np.mean((nodeX,otherXs[0]))
                    midY1 = np.mean((nodeY,otherYs[0]))
                
                    midX2 = np.mean((nodeX,otherXs[1]))
                    midY2 = np.mean((nodeY,otherYs[1]))
                    
                    # get mid points of the sub element
                    midXs = [midX1,midX2]
                    midYs = [midY1,midY2]
                    
                    # now split that sub element into 2 sub sub elements to calcuate
                    # using triangles
                    numSubElements = 2
                
                    # preallocate arrays for use in subsubarea and subarea calculations
                    
                    # for wet area fraction phi
                    # wetDrySubElementList = np.zeros(num_SfcElevs)
                    wetDrySubElementList = cp.zeros(num_SfcElevs)
                    # for grid total water depth H_G
                    # totWatSubElementList = np.zeros(num_SfcElevs)
                    totWatSubElementList = cp.zeros(num_SfcElevs)
                    # for coefficient of friction cf
                    # cfSubElementList = np.zeros(num_SfcElevs)
                    cfSubElementList = cp.zeros(num_SfcElevs)
                    # for wet total water depth H_W
                    # wetTotWatDepthSumList = np.zeros(num_SfcElevs)
                
                    if level0andLevel1:
                        # for Level 1 calculation
                        # rvBottomTermList = np.zeros(num_SfcElevs)
                        rvBottomTermList = cp.zeros(num_SfcElevs)
                        # cadvMiddleTermList = np.zeros(num_SfcElevs)
                        cadvMiddleTermList = cp.zeros(num_SfcElevs)
        
                    countInSubElement = 0
                    
                    # create a list of max elevation to keep track of elevation within 
                    # each subelement then at the end of all the loops take the max
                    # elevation and put it in the maxElevation array
                
                    maxElevationTemp = []
                
                    # loop through each subsubelement then add quantities together
                    # to get averages over entire subelement
                
                    for k in range(numSubElements):
                        
                        # subElementLoopStartTime = time.perf_counter()
                    
                        # get the bounds of the triangle for each subsubArea
                    
                        xCurrSubElement = [nodeX, currEleXCentroid, midXs[k]]
                        yCurrSubElement = [nodeY, currEleYCenteroid, midYs[k]]
                        
                        xySubArray = np.array((xCurrSubElement,yCurrSubElement)).T
                        
                        # convert to meters for area calculations
                        
                        xyCurrSubElementMeters = subgridCalculatormain.projectMeshToMercator(xySubArray[:,1],
                                                                          xySubArray[:,0])
                    
                        
                        # calculate the area of the subsubelement
                        
                        subsubEleArea = subgridCalculatormain.triarea(xyCurrSubElementMeters[0][0],
                                                                      xyCurrSubElementMeters[1][0],
                                                                      xyCurrSubElementMeters[0][1],
                                                                      xyCurrSubElementMeters[1][1],
                                                                      xyCurrSubElementMeters[0][2],
                                                                      xyCurrSubElementMeters[1][2])
                    
                    
                        # add area subelement area calculation to area array
                    
                        area[ele,j] += subsubEleArea
                    
                        # cut down DEM further to save in computation expense
                        
                        # max and min X and Y for the subsub element
                        maxX = (np.max(xCurrSubElement) + 3*xDEMRes) # add a 3 cell buffer
                        minX = (np.min(xCurrSubElement) - 3*xDEMRes)
                    
                        maxY = (np.max(yCurrSubElement) + 3*yDEMRes)
                        minY = (np.min(yCurrSubElement) - 3*yDEMRes)
                    
        
                      # finds rows and column to cut down arrays
                    
                        minRow = np.where(yDEMCoords <= maxY)[0][0]
                        maxRow = np.where(yDEMCoords >= minY)[0][-1]
            
                        minCol = np.where(xDEMCoords >= minX)[0][0]
                        maxCol = np.where(xDEMCoords <= maxX)[0][-1]
                    
            
                        # cut geotiff matrix down for faster processing
                        
                        # cutDownStartTime = time.perf_counter()
                    
                        xCutGeoTiffMatrix2 = xDEM[minRow:maxRow+1,minCol:maxCol+1]#.astype(np.float32)
                        yCutGeoTiffMatrix2 = yDEM[minRow:maxRow+1,minCol:maxCol+1]#.astype(np.float32)
                        zCutGeoTiffMatrix2 = zDEM[minRow:maxRow+1,minCol:maxCol+1]#.astype(np.float32)
                        nCutGeoTiffMatrix2 = nArray[minRow:maxRow+1,minCol:maxCol+1]#.astype(np.float32)
                    
                        # cutDownEndTime = time.perf_counter()
                        
                        # print('Cut Down time = {} s'.format(cutDownEndTime - cutDownStartTime))
                        
                        # for use in determining what cells lie within the subsubelement
                        difCriteria = (np.ones((len(zCutGeoTiffMatrix2[:,0]),
                                            len(zCutGeoTiffMatrix2[0,:])))*areaDif)
                    
                    
                        # mask to find which cells are within a subsubelement
                        # cellsInsideStartTime = time.perf_counter()
                        mask = subgridCalculatormain.isInside(xCurrSubElement[0],yCurrSubElement[0],
                                                              xCurrSubElement[1],yCurrSubElement[1],
                                                              xCurrSubElement[2],yCurrSubElement[2],
                                                              xCutGeoTiffMatrix2,yCutGeoTiffMatrix2,
                                                              difCriteria)
                        
                        # convert mask to cupy array
                        mask = cp.asarray(mask)
                        zCutGeoTiffMatrix2 = cp.asarray(zCutGeoTiffMatrix2)
                        nCutGeoTiffMatrix2 = cp.asarray(nCutGeoTiffMatrix2)
                        surfaceElevations = cp.asarray(surfaceElevations)
                        
                        zCutGeoTiffMatrix2masked = zCutGeoTiffMatrix2[mask]
                        nCutGeoTiffMatrix2masked = nCutGeoTiffMatrix2[mask]
                        
                        # cellsInsideEndTime = time.perf_counter()
                        # print('Isiside time = {} s'.format(cellsInsideEndTime - cellsInsideStartTime))
                        # count how many cells are within subsubelement
                        # countIn = np.count_nonzero(mask)
                        countIn = cp.count_nonzero(mask)
                        
                        # if there are no cells within the element the DEM is too coarse
                        # you must decrease the DEM resolution in this area
                        if countIn == 0:
                        
                            sys.exit('DEM {0} resolution too coarse!'.format(i))
                    
                        # keep track of this for use later
                    
                        countInSubElement += countIn
                    
        ################ BEGIN VECTORIZED CALCULATIONS ##############################
                        # subsubCalcStartTime = time.perf_counter()
                        # create a 3d surface array for use in calculations
                    
                        # tempSurfaceElevArray = cp.ones((len(zCutGeoTiffMatrix2[:,0]),
                        #                               len(zCutGeoTiffMatrix2[0,:]),
                        #                               num_SfcElevs))*surfaceElevations
                        # tempSurfaceElevArray = np.ones((len(zCutGeoTiffMatrix2masked),
                        #                                    num_SfcElevs))*surfaceElevations
                        tempSurfaceElevArray = cp.ones((len(zCutGeoTiffMatrix2masked),
                                                            num_SfcElevs))*surfaceElevations

                        # create a 3d manning array for use in calculations
                    
                        # tempManningArray = nCutGeoTiffMatrix2masked[:,np.newaxis]
                        tempManningArray = nCutGeoTiffMatrix2masked[:,cp.newaxis]
                        # tempManningArray = nCutGeoTiffMatrix2
                    
                        # subtract the bathymetry (2D Array) array from the surface 
                        # elevations (3D array) to get total water depths over the 
                        # sub element

                        # tempTotWatDepthArray = tempSurfaceElevArray -  zCutGeoTiffMatrix2masked[:,np.newaxis]
                        tempTotWatDepthArray = tempSurfaceElevArray -  zCutGeoTiffMatrix2masked[:,cp.newaxis]
                        
                        # max elevation

                        # maxElevationTemp.append(cp.max(zCutGeoTiffMatrix2[mask]))
                        # maxElevationTemp.append(np.max(zCutGeoTiffMatrix2masked))
                        maxElevationTemp.append(cp.max(zCutGeoTiffMatrix2masked))

                        # only look at cells within the subsubelement
                    
                        # tempTotWatDepthArray = tempTotWatDepthArray[mask]
                        # tempManningArray = tempManningArray[mask]
        
                        # find which of these cells are wet
                    
                        # tempWetDryList = tempTotWatDepthArray > 0.0
                        # add some tiny minimum water depth so we dont have cells with
                        # like 10^-6 depths we will use 1 mm to start
                        tempWetDryList = tempTotWatDepthArray > 0.001
                        
                        # count how many cells are wet
                        # tempCountWet = np.count_nonzero(tempWetDryList,axis=0)
                        tempCountWet = cp.count_nonzero(tempWetDryList,axis=0)

                    ################ CALCULATE WET FRACTION #######################
                        # keep track of how many cells are wet for use in averaging 
                        # and wet fraction calulation
                    
                        wetDrySubElementList += tempCountWet
   
                    ###############################################################
        
                    #################### CALCULATING TOTAL WATER DEPTH ############
                    
                        # 0 out any dry cells
                        tempTotWatDepthWetArray = tempTotWatDepthArray * tempWetDryList
                        
                        # integrate the wet total water depths for use in averaging later
                        # totWatSubElementList += cp.nansum(tempTotWatDepthWetArray,axis=0)
                        # totWatSubElementList += np.sum(tempTotWatDepthWetArray,axis=0)
                        totWatSubElementList += cp.sum(tempTotWatDepthWetArray,axis=0)

                        #################### CALCULATING MANNINGS n ###################
                        # find the mannings for only wet areas then nan the rest for 
                        # use in calculations 
                        
                        tempManningWetArray = tempManningArray * tempWetDryList
                        
                        #JLW: i don't think I need to make these nans because we will make
                        # nans in future steps
                        # tempManningWetArray[tempManningWetArray == 0.0] = np.nan
                                                
                        

                        ########### CALCULATE GRID AVERAGED CF FOR MANNING ############
                        
#                         # calulate now for the subsubelement then sum later for use 
#                         # in other calulations
                        tempcf = 9.81*tempManningWetArray**2/tempTotWatDepthWetArray**(1/3)

                        ###############################################################
                        
                        if level0andLevel1:
                            ############ NOW CALCULATE RV FROM KENNEDY ET AL 2019 #########
            
                            # integrate only now and then calculate full rv later
                            # rvBottomTermList += np.nansum(tempTotWatDepthWetArray**(3/2)*\
                            #     (tempcf)**(-1/2),axis = 0)
                            rvBottomTermList += cp.nansum(tempTotWatDepthWetArray**(3/2)*\
                                (tempcf)**(-1/2),axis = 0)

                            ########### CALCULATE Part of C ADVECTION ####################
                            # cadvStartTime = time.perf_counter()
                            # cadvMiddleTermList += np.nansum(tempTotWatDepthWetArray**2/tempcf
                            #                                 ,axis=0)
                            cadvMiddleTermList += cp.nansum(tempTotWatDepthWetArray**2/tempcf
                                                            ,axis=0)
                            

                        # bottomFricEndTime = time.perf_counter()
                        # print('Bottom friction time {} s'.format(bottomFricEndTime-bottomFricStartTime))
                        # finally sum and add cf for use in grid averaged calculation
                        # tempcf = np.nansum(tempcf,axis=0)
                        tempcf = cp.nansum(tempcf,axis=0)
                        cfSubElementList += tempcf
                        # subsubCalcEndTime = time.perf_counter()
                        # print(subsubCalcEndTime-subsubCalcStartTime)
                    # Okay now we can calculate the averaged variable values for each
                    # of the 3 polygonal sub elements within a real adcirc element
                    # subCalcStartTime = time.perf_counter()
                    # go ahead and average wet and grid total water depth for use in calulations
                    # convert back to numpy 
                    # totWatSubElementList = cp.asnumpy(totWatSubElementList)
                    
                    # cadvMiddleTermList = cp.asnumpy(cadvMiddleTermList)
                    # 
                    
                    # countInSubElement = cp.asnumpy(countInSubElement)
                    
                #     # jlw I think the next calulation will give nan regardless
                #     # wetDrySubElementList[wetDrySubElementList == 0] = np.nan
                    wetAvgTotWatDepth = totWatSubElementList/wetDrySubElementList
                    gridAvgTotWatDepth = totWatSubElementList/countInSubElement
                    wetFractionTemp = wetDrySubElementList/countInSubElement
                    cfTemp = cfSubElementList/countInSubElement
                    # rvTemp = wetAvgTotWatDepth/(rvBottomTermList/wetDrySubElementList)
                    cmfTemp = (wetAvgTotWatDepth)*(wetAvgTotWatDepth/(rvBottomTermList/wetDrySubElementList))**2
                    cadvTemp = (1/(wetAvgTotWatDepth)) \
                            *(cadvMiddleTermList/wetDrySubElementList) \
                            *(wetAvgTotWatDepth/(rvBottomTermList/wetDrySubElementList))**2
                    # wetAvgTotWatDepth = cp.asnumpy(wetAvgTotWatDepth)
                    # gridAvgTotWatDepth = cp.asnumpy(gridAvgTotWatDepth)
                    # wetFractionTemp = cp.asnumpy(wetFractionTemp)
                    # cfTemp = cp.asnumpy(cfTemp)
                    # rvBottomTermList = cp.asnumpy(rvBottomTermList)
                    # wetDrySubElementList = cp.asnumpy(wetDrySubElementList)
                    # cfSubElementList = cp.asnumpy(cfSubElementList)
                    # cadvMiddleTermList = cp.asnumpy(cadvMiddleTermList)
   
                    # give grid averaged total water depth
                    # startConvertTime = time.perf_counter()
                    # totWatDepth[ele,j,:] = cp.asnumpy(gridAvgTotWatDepth)
                    totWatDepth[ele,j,:] = cp.ndarray.get(gridAvgTotWatDepth)
                    # endConvertTime = time.perf_counter()
                    # print(endConvertTime-startConvertTime)
                    # give wet averaged total water depth
                    # wetTotWatDepth[ele,j,:] = cp.asnumpy(wetAvgTotWatDepth)
                    wetTotWatDepth[ele,j,:] = cp.ndarray.get(wetAvgTotWatDepth)
                    # get wet area fraction for the subelement 
                    # wetFraction[ele,j,:] = wetDrySubElementList/countInSubElement
                    # wetFraction[ele,j,:] = cp.asnumpy(wetFractionTemp)
                    wetFraction[ele,j,:] = cp.ndarray.get(wetFractionTemp)
                    # get grid averaged coefficient of friction for the subelement
                    ######### JUST TESTING THIS #################
                    # cf[i,j,:] = cfSubElementList/countInSubElement
                    
                    # actually on inspection of implementation, I think I need to solve
                    # for wet averaged cf
                    ########### THIS IS WHAT I NORMALLY USE #################
                    # convert back to numpy
                    
                    # cf[ele,j,:] = cfSubElementList/wetDrySubElementList
                    # cf[ele,j,:] = cp.asnumpy(cfTemp)
                    cf[ele,j,:] = cp.ndarray.get(cfTemp)
                    # get the max elevation of the subgelement
                    maxElevationSubEle[ele,j] = max(maxElevationTemp)
                    
                    if level0andLevel1:
                        
                        # for this rv we average the bottom and top terms separatly over the whole subelement
                        # rv[ele,j,:] = wetAvgTotWatDepth/(rvBottomTermList/wetDrySubElementList)
                        
                        # get corrected bottom friction for level 1 corrections
                        # cmf[i,j,:] = (gridAvgTotWatDepth)*rv[i,j,:]**2 # this is incorrect from the Kennedy et al 2019 paper
                        # cmf[ele,j,:] = (wetAvgTotWatDepth)*rv[ele,j,:]**2
                        # cmf[ele,j,:] = cp.asnumpy(cmfTemp)
                        cmf[ele,j,:] = cp.ndarray.get(cmfTemp)
        
                        # get advection correction for level 1 corrections
                        # cadv[ele,j,:] = (1/(wetAvgTotWatDepth)) \
                        #     *(cadvMiddleTermList/wetDrySubElementList)*rv[ele,j,:]**2
                        # cadv[ele,j,:] = cp.asnumpy(cadvTemp)
                        cadv[ele,j,:] = cp.ndarray.get(cadvTemp)  
                            
                            
                #### ADD MAX ELEVATION ###
                    
                maxElevationEle[ele] = np.max(maxElevationSubEle[ele,:])
                
                countElementLoop += 1
                stopTime = time.perf_counter()
                print("Finished Element {0} of {1} in DEM {2} took {3}".format(countElementLoop,len(elementList),i,stopTime - startTime))
            
                
        # add bottom limit on cf, cmf. and cadv
        cf[cf<0.0025] = 0.0025
        cf[np.isnan(cf)] = 0.0025
        
        if level0andLevel1:
            cmf[cmf<0.0025] = 0.0025   
            cmf[np.isnan(cmf)] = 0.0025
            cadv[np.isnan(cadv)] = 1.0    
            
        wetFraction[np.isnan(wetFraction)] = 0.0
        wetTotWatDepth[np.isnan(wetTotWatDepth)] = 0.0
        totWatDepth[np.isnan(totWatDepth)] = 0.0
        
        
        # put elemental quantities on vertices
        
        wetFractionVertex = np.zeros((numVert,num_SfcElevs))
        wetTotWatDepthVertex = np.zeros((numVert,num_SfcElevs))
        gridTotWatDepthVertex = np.zeros((numVert,num_SfcElevs))
        vertexArea = np.zeros((numVert,num_SfcElevs))
        cfVertex = np.zeros((numVert,num_SfcElevs))
        
        # add max elevation
        maxElevationVertex = np.ones(numVert)*-99999
    
        
        # now fill non subgrid spaces with -99999
        
        wetFractionVertex[np.where(binaryVertexList == 0),:] = -99999
        wetTotWatDepthVertex[np.where(binaryVertexList == 0),:] = -99999
        gridTotWatDepthVertex[np.where(binaryVertexList == 0),:] = -99999
        vertexArea[np.where(binaryVertexList == 0),:] = -99999
        cfVertex[np.where(binaryVertexList == 0),:] = -99999
        maxElevationVertex[np.where(binaryVertexList == 0)] = -99999

        
        if level0andLevel1:
            cmfVertex = np.zeros((numVert,num_SfcElevs))
            cadvVertex = np.zeros((numVert,num_SfcElevs))
            
            # now fill non subgrid spaces with -99999
            
            cmfVertex[np.where(binaryVertexList == 0),:] = -99999
            cadvVertex[np.where(binaryVertexList == 0),:] = -99999
        
        # loop through elements and sum there vertex quantities    
        for j in range(numEle):
            
            if(binaryElementList[j] == 1):
                
                nm1 = meshConnectivity[j,0]
                nm2 = meshConnectivity[j,1]
                nm3 = meshConnectivity[j,2]
                
                phi1 = wetFraction[j,0,:]
                phi2 = wetFraction[j,1,:]
                phi3 = wetFraction[j,2,:]
                
                HW1 = wetTotWatDepth[j,0,:]
                HW2 = wetTotWatDepth[j,1,:]
                HW3 = wetTotWatDepth[j,2,:]
                
                HG1 = totWatDepth[j,0,:]
                HG2 = totWatDepth[j,1,:]
                HG3 = totWatDepth[j,2,:]
                
                cf1 = cf[j,0,:]
                cf2 = cf[j,1,:]
                cf3 = cf[j,2,:]
                
                # add max elevation
                if(maxElevationVertex[nm1] < maxElevationSubEle[j,0]):
                    
                    maxElevationVertex[nm1] = maxElevationSubEle[j,0]
            
                if(maxElevationVertex[nm2] < maxElevationSubEle[j,1]):
                    
                    maxElevationVertex[nm2] = maxElevationSubEle[j,1]
                    
                if(maxElevationVertex[nm3] < maxElevationSubEle[j,2]):
                    
                    maxElevationVertex[nm3] = maxElevationSubEle[j,2]
                
                vertexArea[nm1,:] += area[j,0]
                vertexArea[nm2,:] += area[j,1]
                vertexArea[nm3,:] += area[j,2]
                
                # sum element quantities on vertex 
                # multiply by area to get area weighted average
                
                wetFractionVertex[nm1,:] += phi1 * area[j,0]
                wetFractionVertex[nm2,:] += phi2 * area[j,1]
                wetFractionVertex[nm3,:] += phi3 * area[j,2]
                
                wetTotWatDepthVertex[nm1,:] += HW1 * area[j,0]
                wetTotWatDepthVertex[nm2,:] += HW2 * area[j,1]
                wetTotWatDepthVertex[nm3,:] += HW3 * area[j,2]
                
                gridTotWatDepthVertex[nm1,:] += HG1 * area[j,0]
                gridTotWatDepthVertex[nm2,:] += HG2 * area[j,1]
                gridTotWatDepthVertex[nm3,:] += HG3 * area[j,2]
                
                cfVertex[nm1,:] += cf1 * area[j,0]
                cfVertex[nm2,:] += cf2 * area[j,1]
                cfVertex[nm3,:] += cf3 * area[j,2]
            
                if level0andLevel1:
                
                    cmf1 = cmf[j,0,:]
                    cmf2 = cmf[j,1,:]
                    cmf3 = cmf[j,2,:]
                    
                    cadv1 = cadv[j,0,:]
                    cadv2 = cadv[j,1,:]
                    cadv3 = cadv[j,2,:]
                
                    cmfVertex[nm1,:] += cmf1 * area[j,0]
                    cmfVertex[nm2,:] += cmf2 * area[j,1]
                    cmfVertex[nm3,:] += cmf3 * area[j,2]
                
                    cadvVertex[nm1,:] += cadv1 * area[j,0]
                    cadvVertex[nm2,:] += cadv2 * area[j,1]
                    cadvVertex[nm3,:] += cadv3 * area[j,2]
        
        # now average all of these by the vertex areas
        wetFractionVertex[np.where(binaryVertexList == 1)] = wetFractionVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        wetTotWatDepthVertex[np.where(binaryVertexList == 1)] = wetTotWatDepthVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        gridTotWatDepthVertex[np.where(binaryVertexList == 1)] = gridTotWatDepthVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        cfVertex[np.where(binaryVertexList == 1)] = cfVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        # cfwetVertex = cfwetVertex/vertexArea
        if level0andLevel1:
            
            cmfVertex[np.where(binaryVertexList == 1)] = cmfVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
            # cmfgridVertex = cmfgridVertex/vertexArea
            cadvVertex[np.where(binaryVertexList == 1)] = cadvVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
            # cadvgridVertex = cadvgridVertex/vertexArea
        
        # now we need to check to see if there are any vertices in the subgrid
        # but not connected to any elements in the subgrid 
        
        for i in range(numVert):
            
            # for vertices that are in the subgrid but may or may not be
            # connected to any subgrid elements
            if(binaryVertexList[i]==1):
                # find connected elements
                check = np.where(meshConnectivity == i)[0]
                # check to see if any of the elements are in the subgrid
                inSubgrid = any(binaryElementList[check] == 1)
                
                if(inSubgrid != True):
                    
                    # change the vertex to not be in the subgrid
                    binaryVertexList[i] = 0
                    # change the vertex variables to be -99999
                    wetFractionVertex[i,:] = -99999
                    wetTotWatDepthVertex[i,:] = -99999
                    gridTotWatDepthVertex[i,:] = -99999
                    cfVertex[i,:] = -99999
                    if(level0andLevel1):
                        cmfVertex[i,:] = -99999
                        cadvVertex[i,:] = -99999
                        
            # otherwise if vertex is outside subgrid area
            # just make everything -99999
            else:
                
                # change the vertex variables to be -99999
                wetFractionVertex[i,:] = -99999
                wetTotWatDepthVertex[i,:] = -99999
                gridTotWatDepthVertex[i,:] = -99999
                cfVertex[i,:] = -99999
                if(level0andLevel1):
                    cmfVertex[i,:] = -99999
                    cadvVertex[i,:] = -99999

        # write netcdf file
        
        # if level0andLevel1:
        #     outputFilename = 'E:/Fall2021/varible_dryPhi_tester/test_mesh_1/inputFiles/avgVars_01_increment.nc'
        # else:
        #     outputFilename = 'E:/Fall2021/varible_dryPhi_tester/test_mesh_1/inputFiles/avgVarslevel0only_01_increment.nc'
        
        ncFile = Dataset(outputFilename, mode = 'w', format = 'NETCDF4')
        
        # create dimensions
        
        ncFile.createDimension('numEle',numEle) # element dimension
        ncFile.createDimension('numVert',3) # number of vertices per element
        ncFile.createDimension('numSfcElevs',len(surfaceElevations)) # number of surface elevations
        ncFile.createDimension('numNode',numVert) # number of nodes in mesh
        
        # create variables
        # write variable dimensions transposed because FORTRAN will read them that way
        # only need to do it for the 3D arrays, handle 2 and 1D a different way.
        
        # vertex wet area fraction
        wetFractionVarVertex = ncFile.createVariable('wetFractionVertex',np.float32,
                                            ('numNode','numSfcElevs'))
        # elemental wet area fraction
        wetFractionVarElement = ncFile.createVariable('wetFractionElement',np.float32,
                                                ('numSfcElevs','numVert','numEle'))
        # elemental areas
        areaVar = ncFile.createVariable('area',np.float32,('numEle','numVert'))
        
        # elemental grid averaged total water depth
        totWatDepthVar = ncFile.createVariable('totWatDepth',np.float32,
                                                ('numSfcElevs','numVert','numEle'))
        
        # # vertex minimum wet water depths
        # phiMinDepthVar = ncFile.createVariable('phiMinDepth',np.float32,'numNode')
        
        # surface elevation array
        surfaceElevationsVar = ncFile.createVariable('surfaceElevations',np.float32,'numSfcElevs')
        
        # vertex wet total water depth
        wetTotWatDepthVarVertex = ncFile.createVariable('wetTotWatDepthVertex',np.float32,
                                                          ('numNode','numSfcElevs'))
        
        # vertex grid total water depth
        gridTotWatDepthVarVertex = ncFile.createVariable('gridTotWatDepthVertex',np.float32,
                                                          ('numNode','numSfcElevs'))
        
        # vertex coefficient of friction level 0
        cfVarVertex = ncFile.createVariable('cfVertex',np.float32,
                                            ('numNode','numSfcElevs'))
        
        # variables showing which elements and vertices are contained within
        # the subgrid area
        
        binaryElementListVariable = ncFile.createVariable('binaryElementList',np.int,
                                                          ('numEle'))
        binaryVertexListVariable = ncFile.createVariable('binaryVertexList',np.int,
                                                          ('numNode'))
        
        # write max Elevation
        
        maxElevationEleVariable = ncFile.createVariable('maxElevationElement',np.float32,
                                                      ('numEle'))
        maxElevationVertexVariable = ncFile.createVariable('maxElevationVertex',
                                                            np.float32,
                                                            ('numNode'))
        
        if level0andLevel1:
            # vertex coefficient of friction level 1
            cmfVarVertex = ncFile.createVariable('cmfVertex',np.float32,
                                                  ('numNode','numSfcElevs'))
            
            # elemental advection correction
            cadvVar = ncFile.createVariable('cadv',np.float32,
                                            ('numSfcElevs','numVert','numEle'))
        
        
        wetFractionVarVertex[:,:] = wetFractionVertex
        wetFractionVarElement[:,:,:] = wetFraction.T
        areaVar[:,:] = area
        totWatDepthVar[:,:,:] = totWatDepth.T
        # phiMinDepthVar[:] = phiMinDepth
        surfaceElevationsVar[:] = cp.ndarray.get(surfaceElevations)
        wetTotWatDepthVarVertex[:,:] = wetTotWatDepthVertex
        gridTotWatDepthVarVertex[:,:] = gridTotWatDepthVertex
        cfVarVertex[:,:] = cfVertex
        binaryElementListVariable[:] = binaryElementList
        binaryVertexListVariable[:] = binaryVertexList
        # add max elevation cal
        maxElevationEleVariable[:] = maxElevationEle
        maxElevationVertexVariable[:] = maxElevationVertex
        # cfVarElement[:,:,:] = cf
        if level0andLevel1:
            
            cmfVarVertex[:,:] = cmfVertex
            # cmfVarElement[:,:,:] = cmf
            cadvVar[:,:,:] = cadv.T
        
        ncFile.close()
        
        endTot = time.perf_counter()
        
        print('Total Time = {} minutes'.format((endTot-startTot)/60))
        
        # return phi1011Check, phi1057Check, phi1056Check, phi1013Check, phi1010Check, phi1053Check
        
        # # try using a patch plot
         
        # plot vertex quantities
        
        # fig,ax = plt.subplots(3,1,figsize=(12,9))
        
        # ax[0].plot(surfaceElevations,cfVertex[186,:],label='<Cf>W')
        # # ax[0].plot(surfaceElevations,cfwetVertex[490,:],label='<Cf>W')
        # if level0andLevel1:
        #     ax[0].plot(surfaceElevations,cmfVertex[186,:],label='<H>W*Rv^2')
        #     # ax[0].plot(surfaceElevations,cmfgridVertex[490,:],label='<H>G*Rv^2')
        #     ax[0].set_title('cf vs cmf')
        # else:
        #     ax[0].set_title('cf')
            
        # ax[0].legend()
        
        # ax[1].plot(surfaceElevations,wetFractionVertex[186,:])     
        # ax[1].set_title('phi')   
        
        # if level0andLevel1:
        #     ax[2].plot(surfaceElevations,cadvVertex[186,:],label='wetavg')   
        #     # ax[2].plot(surfaceElevations,cadvgridVertex[490,:],label='gridavg middleG')   
        #     ax[2].set_title('cadv')  
        #     ax[2].set_ylim((0,2))
            
        # cmap = cmocean.cm.rain
        # colorLevels = np.arange(0,1.05,0.05)
        # cbarTicks = np.arange(0,1.1,0.1)
            
        # adcircToolsModule.resultsPlotter.timestep63plotter(wetFractionVertex[:,200],mesh[1],cmap,
        #                               colorLevels,cbarTicks,'phi Vertex')           
                    
