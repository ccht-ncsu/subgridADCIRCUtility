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
        
    def readManning(manningsnFilename):
        
        # read mannings n file and store values in a dictionary
        
        manningsValues =  open(manningsnFilename,'r')
        
        manningsnTable = {}
        
        for value in manningsValues:
            
            line = value.split()
            
            manningsnTable[int(line[0])] = float(line[1])
            
        return manningsnTable
    
######################## GET TEXT FILE WITH MESH RESOLUTION ###################
    
    def meshResolution(meshObject,outputFilename):
        
        import numpy  as np
        
        # get the mesh connectivity 
        
        meshConnectivity = meshObject[1].triangles
        
        # get the x and y of mesh Vertices
        
        allVertLon = meshObject[0]['Longitude']
        allVertLat = meshObject[0]['Latitude']
        
        # create empty list of lists
        
        connectedVertices = [ [] for _ in range(meshObject[2]) ]
        
        # create array to hold distances
        distToConnectedVertices = np.zeros(meshObject[2])
       
        # get vertex connectivity
        
        for i in range(meshObject[3]):
            
            currMeshConn = meshConnectivity[i]
            
            nm0 = currMeshConn[0]
            nm1 = currMeshConn[1]
            nm2 = currMeshConn[2]
            
            connectedVertices[nm0] = connectedVertices[nm0] + [currMeshConn[1],currMeshConn[2]]
            connectedVertices[nm1] = connectedVertices[nm1] + [currMeshConn[0],currMeshConn[2]]
            connectedVertices[nm2] = connectedVertices[nm2] + [currMeshConn[0],currMeshConn[1]]
        
        
        # get unique values and calulate distances
        
        for i in range(meshObject[2]):

            # get unique values 
            connectedVertices[i] = np.array(connectedVertices[i])
            connectedVertices[i] = np.unique(connectedVertices[i])
            
            # get x y of vertex of interest
            
            vertLon = allVertLon[i]
            vertLat = allVertLat[i]
            
            # now get the x y of the connected vertices
            
            conVertLon = allVertLon[connectedVertices[i]]
            conVertLat = allVertLat[connectedVertices[i]]
        
            # now calculate distances
            
            conVertDistances = np.sqrt((conVertLon - vertLon)**2 +
                                       (conVertLat - vertLat)**2) * 111/0.001 # convert to meters
            
            # average these
            
            convertDistancesMean = np.mean(conVertDistances)
            
            # now add this to larger array
        
            distToConnectedVertices[i] = convertDistancesMean
        
        with open(outputFilename,'w+') as meshRes:
            
            meshRes.write('Averaged distance surrounding each vertex\n')
            
            for i in range(meshObject[2]):
                
                meshRes.write(str(i)+'\t'+str(distToConnectedVertices[i])+'\n')        



########################## SUBGRID CALULATOR OLD SCHEME ################################

    def calculateSubgridCorrectionOld(controlFilename):
        
        import sys
        import netCDF4 as nc
        import numpy as np
        import time
        import matplotlib.pyplot as plt
        import re
 
        # first read in the control file
        
        controlFile = controlFilename
        
        with open(controlFile) as ctrF:
            
            ctrF.readline()
            # change to shintaros r.strip with re to allow for spaces
            line = ctrF.readline().rstrip()
            line = re.split(' *= *',line)
            # get output file name
            outputFilename = line[1]
            line = ctrF.readline().rstrip()
            line = re.split(' *= *',line)
            # get mesh filename
            meshFilename = line[1]
            # read in mannings stuff
            line = ctrF.readline().rstrip()
            line = re.split(' *= *',line)
            defaultManning = line[1] # if true just use the manning table in the code
            if defaultManning == 'False': # otherwise we need to read a table in
                line = ctrF.readline().rstrip()
                line = re.split(' *= *',line)
                manningsnTableFilename = line[1] # get the mannings n table filename
            # now read in the elevation array
            line = ctrF.readline().rstrip()
            line = re.split(' *= *',line)
            minSurElev = float(line[1])
            line = ctrF.readline().rstrip()
            line = re.split(' *= *',line)
            maxSurElev = float(line[1])
            line = ctrF.readline().rstrip()
            line = re.split(' *= *',line)
            elevDisc = float(line[1])
            line = ctrF.readline().rstrip()
            line = re.split(' *= *',line)
            numDEMs = int(line[1])
            # get list of elevation datasets
            demFilenameList = []
            for i in range(numDEMs):
                line = ctrF.readline().rstrip()
                line = re.split(' *= *',line)
                demFilenameList.append(line[0])
            line = ctrF.readline().rstrip()
            line = re.split(' *= *',line)
            numLCs = int(line[1])
            # get list of landcover datasets
            landcoverFilenameList = []
            for i in range(numLCs):
                line = ctrF.readline().rstrip()
                line = re.split(' *= *',line)
                landcoverFilenameList.append(line[0])
        
        # first lets import the mesh 
        
        mesh = subgridCalculatormain.readMesh(meshFilename)
        
        meshConnectivity = mesh[1].triangles
        meshLon = np.asarray(mesh[0]['Longitude']).astype('float64')
        meshLat = np.asarray(mesh[0]['Latitude']).astype('float64')
        numNode = mesh[2]
        numEle = mesh[3]
        
                
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
            elevationData = subgridCalculatormain.importDEMv2(demFilenameList[i])
            # get bathy/topo elevations
            zDEMTemp = elevationData[0]
            # resolution in the x direction
            xDEMResTemp = elevationData[1]
            # resolution in the y direction
            yDEMResTemp = -1*elevationData[2]
            # x coordinates of DEM in 1D array
            xDEMCoordsTemp = elevationData[3]
            # y coordinates of DEM in 1D array
            yDEMCoordsTemp = elevationData[4]
            # deallocate DEM data 
            elevationData = None
            
            # get dem bounds to determine what dem to use when looping
            
            elevationDict["bounds%s"%i] = [np.min(xDEMCoordsTemp),
                                            np.max(xDEMCoordsTemp),
                                            np.min(yDEMCoordsTemp),
                                            np.max(yDEMCoordsTemp)]
            
            print('Finished reading DEM {0}.'.format(i))
        
        # deallocate arrays
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
        
        # deallocate Xs and Ys
        
        Xs = None
        Ys = None
        
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
        
        for i in range(len(demFilenameList)):
        
            # fine which elements are in the DEM
            totalEleInfoTable[:,5] = ((totalEleInfoTable[:,1]>elevationDict["bounds%s"%i][0])
                                      & ((totalEleInfoTable[:,2])<elevationDict["bounds%s"%i][1])
                                      & ((totalEleInfoTable[:,3])>elevationDict["bounds%s"%i][2])
                                      & ((totalEleInfoTable[:,4])<elevationDict["bounds%s"%i][3]))
        
            whichAreInside = list(np.where(totalEleInfoTable[:,5] == 1)[0])
            elementDict["DEM%s"%i] = totalEleInfoTable[whichAreInside,0]
            
            # create a list of elements within subgrid area
            
            containedElementList0Index.append(whichAreInside)
            
            # create list of vertices within subgrid area
            
            totalVertInfoTable[:,3] = ((totalVertInfoTable[:,1]>elevationDict["bounds%s"%i][0])
                                      & ((totalVertInfoTable[:,1])<elevationDict["bounds%s"%i][1])
                                      & ((totalVertInfoTable[:,2])>elevationDict["bounds%s"%i][2])
                                      & ((totalVertInfoTable[:,2])<elevationDict["bounds%s"%i][3]))
            
            whichAreInside = list(np.where(totalVertInfoTable[:,3]==1)[0])
            
            containedVertexList0Index.append(whichAreInside)
            
        
        # make sure each dem only has unique element numbers
        for i in range(len(elementDict)):
            
            currContainedElements = elementDict['DEM%s'%i]
            
            for j in range(i+1,len(elementDict)):
                
                currOtherContainedElements = elementDict['DEM%s'%j]
                overLappingElement = np.intersect1d(currContainedElements,
                                                    currOtherContainedElements,
                                                    return_indices=True)
                # now delete indices in the other dems that correspond with 
                # ones that have already been read
                elementDict['DEM%s'%j] = np.delete(elementDict['DEM%s'%j],
                                                   overLappingElement[2],
                                                   axis=0)

        
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
        
        # surface elevation array for calcuations
        surfaceElevations = np.round(np.arange(minSurElev,maxSurElev+elevDisc,
                                               elevDisc),2).astype('float32') 
        
        # preallocate necessary arrays
        numEle = mesh[3]
        nodesPerEle = 3
        num_SfcElevs = len(surfaceElevations)
        
        wetFraction = np.zeros((numEle,nodesPerEle,num_SfcElevs)).astype(np.float32)
        
        area = np.zeros((numEle,nodesPerEle)).astype(np.float32)
        
        totWatDepth = np.zeros((numEle,nodesPerEle,num_SfcElevs)).astype(np.float32)
        
        wetTotWatDepth = np.zeros((numEle,nodesPerEle,num_SfcElevs)).astype(np.float32)
        
        cf = np.zeros((numEle,nodesPerEle,num_SfcElevs)).astype(np.float32)
        
        rv = np.zeros((numEle,nodesPerEle,num_SfcElevs)).astype(np.float32)
                    
        cmf = np.zeros((numEle,nodesPerEle,num_SfcElevs)).astype(np.float32)
                    
        cadv = np.zeros((numEle,nodesPerEle,num_SfcElevs)).astype(np.float32)
        
        # do not use this in calculations anymore JLW 01302024
        # maxElevationEle = np.zeros(numEle).astype(np.float32) # find highest elevation in each element for use in variable phi
        
        # maxElevationSubEle = np.zeros((numEle,3)).astype(np.float32) # find the highest elevation in each subElement
        
        # now fill the rows and columns of non subgrid vertices and elements with -99999
        
        wetFraction[np.where(binaryElementList == 0),:,:] = -99999
        area[np.where(binaryElementList == 0)] = -99999
        totWatDepth[np.where(binaryElementList == 0),:,:] = -99999
        wetTotWatDepth[np.where(binaryElementList == 0),:,:] = -99999
        cf[np.where(binaryElementList == 0),:,:] = -99999
        
        # do not use this in calculations anymore JLW 01302024
        # # fill max elevation
        # maxElevationEle[np.where(binaryElementList == 0)] = -99999
        # maxElevationSubEle[np.where(binaryElementList == 0)] = -99999
        
        # now lets loop through the elements and do some calculations!
        
        # ignore true divide error from wet averaged quantities
        np.seterr(invalid='ignore')
        
        startElemental = time.time()
        
        for i in range(len(demFilenameList)):
            
            # reading in DEM again
            # all variables the same as before
            elevationData = subgridCalculatormain.importDEMv2(demFilenameList[i])
            landcoverData = subgridCalculatormain.importDEMv2(landcoverFilenameList[i])
            
            bathyTopo = elevationData[0].astype('float32')
            lonRes = elevationData[1]
            latRes = -1*elevationData[2]
            lon = elevationData[3]
            lat = elevationData[4]
            elevationData = None # deallocate 
            manningsn = landcoverData[0].astype('float32') # array of mannings n values
            landcoverData = None # deallocate
            
            if defaultManning:
            
                landCoverToManning = {0:0.02, 2:0.15, 3:0.10, 4:0.05, 5:0.02,
                                      6:0.037, 7:0.033, 8:0.034, 9:0.1, 10:0.11,
                                      11:0.1, 12:0.05, 13:0.1, 14:0.048, 15:0.045,
                                      16:0.1, 17:0.048, 18:0.045, 19:0.04,
                                      20:0.09, 21:0.02, 22:0.015, 23:0.015, 
                                      24:0.09, 25:0.01}
                
            else:
                
                landCoverToManning = subgridCalculatormain.readManning(manningsnTableFilename)

            landCoverValues = [0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25]
            # covert values
            
            for value in landCoverValues:
            
                manningsn[manningsn == value] = landCoverToManning[value]
            
            # set nan values to 0.02
            
            manningsn[np.isnan(manningsn)] = 0.02
            
            # get a list of elements within this DEM
            elementList = elementDict["DEM%s"%i]
            
            countElementLoop = 0
        
            for ele in elementList:
            # for ele in range(2):
                
                start = time.time()
                
                # cast ele to integer
                ele = int(ele)
                
                # get vertex numbers
            
                nm0 = meshConnectivity[ele,0]
                nm1 = meshConnectivity[ele,1]
                nm2 = meshConnectivity[ele,2]
                
                # put this in a list
                
                nodeNumbers = [nm0,nm1,nm2]
                
                # now get vertex coordinates
                
                nodeLon = meshLon[nodeNumbers]
                nodeLat = meshLat[nodeNumbers]
            
                # now get the centroid of the element
                
                centroidLon = np.mean(nodeLon)
                centroidLat = np.mean(nodeLat)
                
                # now get the mid points 
                
                midLon = []
                midLat = []
                
                midIndex = [[0,1],[1,2],[2,0]]

                for j in range(3):
                    
                    midLon.append(np.mean(nodeLon[midIndex[j]]))
                    midLat.append(np.mean(nodeLat[midIndex[j]]))
                
                
                subAreaPerimeter = np.zeros((2,4))
                
                perIndex = [[0,0,2],
                            [0,1,1],
                            [1,2,2]]
            
                for j in range(3):
                    
                    subAreaPerimeter[:,:] = ([[centroidLon,midLon[perIndex[j][2]],
                                                  nodeLon[perIndex[j][1]],
                                                  midLon[perIndex[j][0]]],
                                              [centroidLat,midLat[perIndex[j][2]],
                                                nodeLat[perIndex[j][1]],
                                                midLat[perIndex[j][0]]]])

                    
                    # cut down dem and landcover to each sub area
                    
                    minEleLon = np.min(subAreaPerimeter[0,:])
                    minEleLat = np.min(subAreaPerimeter[1,:])
                    maxEleLon = np.max(subAreaPerimeter[0,:])
                    maxEleLat = np.max(subAreaPerimeter[1,:])
                    
                    rows = list(np.where((lon > minEleLon) * (lon < maxEleLon))[0])
                    minRow = np.min(rows)
                    maxRow = np.max(rows)
                    demLonCut = lon[rows]
                    cols = list(np.where((lat > minEleLat) * (lat < maxEleLat))[0])
                    minCol = np.min(cols)
                    maxCol = np.max(cols)
                    demLatCut = lat[cols]
                    demBathyTopoCut = bathyTopo[minCol:maxCol+1,minRow:maxRow+1]
                    manningsnCut = manningsn[minCol:maxCol+1,minRow:maxRow+1]
                    
                    lonGrid,latGrid = np.meshgrid(demLonCut,demLatCut)
                    
                    # split into 2 triangles
                
                    tri0 = subAreaPerimeter[:,:3]
                    
                    # convert to meters
                    tri0Meters = subgridCalculatormain.projectMeshToMercator(tri0[1,:],
                                                                      tri0[0,:])
                    # print(tri0)
                    tri0Area = subgridCalculatormain.triarea(tri0Meters[0][0], 
                                                                tri0Meters[1][0],
                                                                tri0Meters[0][1], 
                                                                tri0Meters[1][1],
                                                                tri0Meters[0][2], 
                                                                tri0Meters[1][2])

                    
                    insideTri0 = subgridCalculatormain.isInside(tri0[0,0], tri0[1,0],
                                                                    tri0[0,1], tri0[1,1],
                                                                    tri0[0,2], tri0[1,2],
                                                                    lonGrid, latGrid, 0.00000001)
                    
                    tri1 = subAreaPerimeter[:,[0,2,3]]
                    # print(tri1)
                    tri1Meters = subgridCalculatormain.projectMeshToMercator(tri1[1,:],
                                                                      tri1[0,:])

                    
                    tri1Area = subgridCalculatormain.triarea(tri1Meters[0][0], 
                                                                tri1Meters[1][0],
                                                                tri1Meters[0][1], 
                                                                tri1Meters[1][1],
                                                                tri1Meters[0][2], 
                                                                tri1Meters[1][2])
                    
                    
                    insideTri1 = subgridCalculatormain.isInside(tri1[0,0], tri1[1,0],
                                                                    tri1[0,1], tri1[1,1],
                                                                    tri1[0,2], tri1[1,2],
                                                                    lonGrid, latGrid, 0.00000001)
                    
                    # now combine the two triangles and find the points inside
                    
                    insideSubElement = np.logical_or(insideTri0,insideTri1)
                    
                    # count the number of subgrid cells within the subelement
                    
                    cellsInSubElement = np.count_nonzero(insideSubElement)
                    
                    # if there are no cells within the element the DEM is too coarse
                    # you must decrease the DEM resolution in this area
                    if cellsInSubElement == 0:
    
                        sys.exit('DEM {0} resolution too coarse!'.format(i))
                    
                    # get just he bathy topo inside the sub element 
                    
                    bathyTopoInsideSubElement = demBathyTopoCut*insideSubElement
                    
                    # set 0 values to nan for calculations
                    
                    bathyTopoInsideSubElement[bathyTopoInsideSubElement==0] = np.nan

                    # do not use this in calculations anymore JLW 01302024
                    # get maximum elevation inside the sub element
                    # 
                    # maxElevationSubEle[ele,j] = np.nanmax(bathyTopoInsideSubElement)
                    
                    # get area of sub element
                    
                    area[ele,j] = tri0Area + tri1Area
                    
                    # print(area[ele,j])
                    
                    # vectorize for surface elevation loop to speed up calcs
                    # remove cells not inside sub element which will flatten the array
                    bathyTopoInsideSubElementNoNaN = bathyTopoInsideSubElement[~np.isnan(bathyTopoInsideSubElement)]
                    manningsnCutNoNaN = manningsnCut[~np.isnan(bathyTopoInsideSubElement)]
                    
                    # get the total water depth at each surface elevation
                    
                    temptotWatDepth =  surfaceElevations[:,None] - bathyTopoInsideSubElementNoNaN
                    
                    # count the number of wet cells
                    
                    wetCellsInSubArea = temptotWatDepth > 0.0001
                    # wetCellsInSubArea = temptotWatDepth > 0.01 # changed to 1 cm
                    
                    wetCellsInSubAreaCount = np.sum(wetCellsInSubArea,axis=1)

                    # now set tot water depth of dry cells to nan
                    
                    temptotWatDepth[temptotWatDepth < 0.0001] = np.nan
                    # temptotWatDepth[temptotWatDepth < 0.01] = np.nan # changed to 1 cm
                    
                    # add to wet frac array
                    
                    wetFraction[ele,j,:] = wetCellsInSubAreaCount/cellsInSubElement
                    
                    # add to total water depth array
                                    
                    totWatDepth[ele,j,:] = np.nansum(temptotWatDepth,axis=1)/cellsInSubElement
                    
                    # get wet total water depth and coefficient of friction
                    
                    # find the mannings for only wet areas then 0 the rest for 
                    # use in calculations 
                    
                    manningsnCutNoNaNWet = manningsnCutNoNaN * wetCellsInSubArea
                    
                    tempcf = (9.81*manningsnCutNoNaNWet**2)/(temptotWatDepth**(1/3))
                    
                    # get wet averaged total water depth
                    
                    wetTotWatDepth[ele,j,:] = np.nansum(temptotWatDepth,axis=1)/wetCellsInSubAreaCount
                    
                    # get bottom friction
                    cf[ele,j,:] = np.nansum(tempcf,axis=1)/wetCellsInSubAreaCount # this is correct I need <Cf>W
                    # cf[ele,j,:] = np.nansum(tempcf,axis=1)/cellsInSubElement
                    
                    # get rv for advection correction and bottom friction correction
                    
                    rv[ele,j,:] = wetTotWatDepth[ele,j,:]/(np.nansum((temptotWatDepth**(3/2))*(tempcf**(-1/2)),axis=1)/wetCellsInSubAreaCount)
                    
                    # get advection correction
                    
                    cadv[ele,j,:] = (1/wetTotWatDepth[ele,j,:])*(np.nansum(temptotWatDepth**2/tempcf,axis=1)/wetCellsInSubAreaCount)*rv[ele,j,:]**2
                    
                    # get corrected bottom friction for level 1 corrections
                    
                    cmf[ele,j,:] = wetTotWatDepth[ele,j,:]*rv[ele,j,:]**2 # this is correct I need <H>W * Rv**2
                # do not use this in calculations anymore JLW 01302024
                # maxElevationEle[ele] = np.max(maxElevationSubEle[ele,:])
            
                end = time.time()
                
                countElementLoop += 1
             
                # now get the quadrilaterals by each node
            
                print("Finished Element {0} of {1} in DEM {2} took {3}".format(countElementLoop,len(elementList),i,end - start))
            
        endElemental = time.time()
        
        print("Elemental Calculations Took {}s".format(endElemental - startElemental))
                        
        cf[cf<0.0025] = 0.0025
        cmf[cmf<0.0025] = 0.0025   
        cf[np.isnan(cf)] = 0.0025
        cmf[np.isnan(cmf)] = 0.0025
        cadv[np.isnan(cadv)] = 1.0  
        
        wetFraction[np.isnan(wetFraction)] = 0.0
        wetTotWatDepth[np.isnan(wetTotWatDepth)] = 0.0
        totWatDepth[np.isnan(totWatDepth)] = 0.0
        
        # return wetTotWatDepth
        
        # now put sub-elemental quantities on the vertices
        
        wetFractionVertex = np.zeros((numNode,num_SfcElevs))
        wetTotWatDepthVertex = np.zeros((numNode,num_SfcElevs))
        gridTotWatDepthVertex = np.zeros((numNode,num_SfcElevs))
        vertexArea = np.zeros((numNode,num_SfcElevs))
        cfVertex = np.zeros((numNode,num_SfcElevs))
        cmfVertex = np.zeros((numNode,num_SfcElevs))
        # cadvVertex = np.zeros((numNode,num_SfcElevs))
        # do not use this in calculations anymore JLW 01302024
        # maxElevationVertex = np.ones(numNode)*-99999
        
        # now fill non subgrid spaces with -99999
        
        wetFractionVertex[np.where(binaryVertexList == 0),:] = -99999
        wetTotWatDepthVertex[np.where(binaryVertexList == 0),:] = -99999
        gridTotWatDepthVertex[np.where(binaryVertexList == 0),:] = -99999
        vertexArea[np.where(binaryVertexList == 0),:] = -99999
        cfVertex[np.where(binaryVertexList == 0),:] = -99999
        # do not use this in calculations anymore JLW 01302024
        # maxElevationVertex[np.where(binaryVertexList == 0)] = -99999
        
        # loop through the elements and sum sub-element quantities surrounding each vertex
        start = time.time()
        
        for i in range(numEle):
            
            if(binaryElementList[i] == 1):
            
                nm0 = meshConnectivity[i,0]
                nm1 = meshConnectivity[i,1]
                nm2 = meshConnectivity[i,2]
                
                phi0 = wetFraction[i,0,:]
                phi1 = wetFraction[i,1,:]
                phi2 = wetFraction[i,2,:]
                
                HW0 = wetTotWatDepth[i,0,:]
                HW1= wetTotWatDepth[i,1,:]
                HW2 = wetTotWatDepth[i,2,:]
                
                HG0 = totWatDepth[i,0,:]
                HG1 = totWatDepth[i,1,:]
                HG2 = totWatDepth[i,2,:]
                
                cf0 = cf[i,0,:]
                cf1 = cf[i,1,:]
                cf2 = cf[i,2,:]
                
                # do not use this in calculations anymore JLW 01302024
                # # add max elevation
                # if(maxElevationVertex[nm0] < maxElevationSubEle[i,0]):
                    
                #     maxElevationVertex[nm0] = maxElevationSubEle[i,0]
            
                # if(maxElevationVertex[nm1] < maxElevationSubEle[i,1]):
                    
                #     maxElevationVertex[nm1] = maxElevationSubEle[i,1]
                    
                # if(maxElevationVertex[nm2] < maxElevationSubEle[i,2]):
                    
                #     maxElevationVertex[nm2] = maxElevationSubEle[i,2]
                
                vertexArea[nm0,:] += area[i,0]
                vertexArea[nm1,:] += area[i,1]
                vertexArea[nm2,:] += area[i,2]
                
                # sum element quantities on vertex 
                # multiply by area to get area weighted average
                
                wetFractionVertex[nm0,:] += phi0 * area[i,0]
                wetFractionVertex[nm1,:] += phi1 * area[i,1]
                wetFractionVertex[nm2,:] += phi2 * area[i,2]
                
                wetTotWatDepthVertex[nm0,:] += HW0 * area[i,0]
                wetTotWatDepthVertex[nm1,:] += HW1 * area[i,1]
                wetTotWatDepthVertex[nm2,:] += HW2 * area[i,2]
                
                gridTotWatDepthVertex[nm0,:] += HG0 * area[i,0]
                gridTotWatDepthVertex[nm1,:] += HG1 * area[i,1]
                gridTotWatDepthVertex[nm2,:] += HG2 * area[i,2]
                
                cfVertex[nm0,:] += cf0 * area[i,0]
                cfVertex[nm1,:] += cf1 * area[i,1]
                cfVertex[nm2,:] += cf2 * area[i,2]
                
                cmf0 = cmf[i,0,:]
                cmf1 = cmf[i,1,:]
                cmf2 = cmf[i,2,:]
            
                cmfVertex[nm0,:] += cmf0 * area[i,0]
                cmfVertex[nm1,:] += cmf1 * area[i,1]
                cmfVertex[nm2,:] += cmf2 * area[i,2]
        
        # now average all of these by the vertex areas
        wetFractionVertex[np.where(binaryVertexList == 1)] = wetFractionVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        wetTotWatDepthVertex[np.where(binaryVertexList == 1)] = wetTotWatDepthVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        gridTotWatDepthVertex[np.where(binaryVertexList == 1)] = gridTotWatDepthVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        cfVertex[np.where(binaryVertexList == 1)] = cfVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        cmfVertex[np.where(binaryVertexList == 1)] = cmfVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]

        end = time.time()
        
        print('Put elemental quantities to vertices {} s'.format(end-start))

        # now we need to check to see if there are any vertices in the subgrid
        # but not connected to any elements in the subgrid 
        start = time.time()
        
        meshConnectivityInSubgrid = meshConnectivity[containedElementList0Index]
        # now just get the unique vertex numbers to reduce the size of this array
        meshConnectivityInSubgrid = np.unique(meshConnectivityInSubgrid)
        # now get an array of vertices not contained in the subgrid
        vertNotInSubgrid = np.delete(totalVertInfoTable[:,0].astype(int),
                                     meshConnectivityInSubgrid)
        # change the vertex to not be in the subgrid
        binaryVertexList[vertNotInSubgrid] = 0
        # change the vertex variables to be -99999
        wetFractionVertex[vertNotInSubgrid,:] = -99999
        wetTotWatDepthVertex[vertNotInSubgrid,:] = -99999
        gridTotWatDepthVertex[vertNotInSubgrid,:] = -99999
        cfVertex[vertNotInSubgrid,:] = -99999
        cmfVertex[vertNotInSubgrid,:] = -99999
                
        desiredPhiList = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        
        end = time.time()
        print('Checked if vertex in subgrid {} s'.format(end-start))
        # transpose the dimensions of the elemental subgrid arrays for writing 
        # netCDF
        start = time.time()
        
        wetFraction = wetFraction.T
        totWatDepth = totWatDepth.T
        cadv = cadv.T
        
        
        # create empty arrays for the reduced tables
        
        depthsEleForLookup = np.zeros((11,3,len(wetFraction[0,0,:])))
        HEleForLookup = np.zeros((11,3,len(wetFraction[0,0,:])))
        cadvForLookup = np.zeros((11,3,len(wetFraction[0,0,:])))
                            
        minSurfElev = np.min(surfaceElevations)
        maxSurfElev = np.max(surfaceElevations)
        
        end = time.time()
        print('Finished prepping for Reduction took {} s'.format(end-start))

        start = time.time()
        
        for i in range(numEle):
            
            element = i
            
            for j in range(3):
                
                vert = j
                
                currPhiArray = wetFraction[:,vert,element]
            
                # find where currPhiArray is equal to 0
                equalTo0 = np.where(currPhiArray == 0.0)[0]
                
                if(len(equalTo0)!=0): # if 0.0 exists in the array
                
                    depthsEleForLookup[0,vert,element] = surfaceElevations[equalTo0[-1]]
                    HEleForLookup[0,vert,element] = totWatDepth[equalTo0[-1],vert,element]
                    cadvForLookup[0,vert,element] = cadv[equalTo0[-1],vert,element]
                    
                else: # so if it never gets fully dry set everything to the value corresponding to the first surface elevations
                
                    depthsEleForLookup[0,vert,element] = surfaceElevations[0]
                    HEleForLookup[0,vert,element] = totWatDepth[0,vert,element]
                    cadvForLookup[0,vert,element] = cadv[0,vert,element]
                    
                # now check for when phi == 1.0 and find exactly where that is
                
                equalTo1 = np.where(currPhiArray == 1.0)[0]
                
                if(len(equalTo1)!=0): # if 1.0 exists in the array

                    depthsEleForLookup[-1,vert,element] = surfaceElevations[equalTo1[0]]
                    HEleForLookup[-1,vert,element] = totWatDepth[equalTo1[0],vert,element]
                    cadvForLookup[-1,vert,element] = cadv[equalTo1[0],vert,element]
                    
                else: # if there is nothing that is equal to 1 (so never gets fully wet, just set everything to correspind to the last surface elevation)
                
                    depthsEleForLookup[-1,vert,element] = surfaceElevations[-1]
                    HEleForLookup[-1,vert,element] = totWatDepth[-1,vert,element]
                    cadvForLookup[-1,vert,element] = cadv[-1,vert,element]
                    
                # now for everything else
                
                for k in range(1,len(desiredPhiList)-1):
                
                    desiredPhi = desiredPhiList[k]
                    greaterThan = np.where(currPhiArray > desiredPhi)[0]
                    
                    if(len(greaterThan)==0):  # so if nothing in the currPhiArray is greater than the desired phi 
                    
                    # set everything to correspond to the last surface elevation
                    
                        depthsEleForLookup[k,vert,element] = surfaceElevations[-1]
                        HEleForLookup[k,vert,element] = totWatDepth[-1,vert,element]
                        cadvForLookup[k,vert,element] = cadv[-1,vert,element]
                        
                    elif(greaterThan[0] == 0): # so if the first currphi index is greater than the desired phi 
                    
                    # set everything to correspond to the first surfaceelevation

                        depthsEleForLookup[k,vert,element] = surfaceElevations[0]
                        HEleForLookup[k,vert,element] = totWatDepth[0,vert,element]
                        cadvForLookup[k,vert,element] = cadv[0,vert,element]
            
                        
                    else: # this is where we interpolate 
                        
                        greaterThan = greaterThan[0]
                        lessThan = greaterThan - 1
                  
                        depthsEleForLookup[k,vert,element] = (((desiredPhi - currPhiArray[lessThan])
                                                      /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                      *(surfaceElevations[greaterThan] - surfaceElevations[lessThan])
                                                      + (surfaceElevations[lessThan]))
                        
                        HEleForLookup[k,vert,element] = (((desiredPhi - currPhiArray[lessThan])
                                                      /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                      *(totWatDepth[greaterThan,vert,element] - totWatDepth[lessThan,vert,element])
                                                      + (totWatDepth[lessThan,vert,element]))
                        
                        cadvForLookup[k,vert,element] = (((desiredPhi - currPhiArray[lessThan])
                                                      /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                      *(cadv[greaterThan,vert,element] - cadv[lessThan,vert,element])
                                                      + (cadv[lessThan,vert,element]))
                    
                  
        # fill the elements that are not contained in the subgrid region with -99999
        depthsEleForLookup[:,:,np.where(binaryElementList==0)[0]] = -99999
        HEleForLookup[:,:,np.where(binaryElementList==0)[0]] = -99999
        cadvForLookup[:,:,np.where(binaryElementList==0)[0]] = -99999
        # deallocate arrays
        totWatDepth = None
        cadv = None
        wetFraction = None
        end = time.time()
        print('Reduction of partially wet elements took {} s'.format(end-start))
                    
#### NOW WE NEED TO REPEAT THIS FOR THE VERTICES ########
        start = time.time() 
        depthsVertForLookup = np.zeros((len(wetFractionVertex[:]),11))
        HGVertForLookup = np.zeros((len(wetFractionVertex[:]),11))
        HWVertForLookup = np.zeros((len(wetFractionVertex[:]),11))
        cfVertForLookup = np.zeros((len(wetFractionVertex[:]),11))
        cmfVertForLookup = np.zeros((len(wetFractionVertex[:]),11))
            
        for i in range(numNode):
            
            vert = i
            currPhiArray = wetFractionVertex[vert,:]
            
            # make sure that the phi array also gets fully wet and then proceed
            # otherwise just skip
            
            # for phi == 0 you want to find exactly where that is in the currPhiArray
            equalTo0 = np.where(currPhiArray == 0.0)[0]
            
            if(len(equalTo0)!=0): # if 0.0 exists in the array
            
                depthsVertForLookup[vert,0] = surfaceElevations[equalTo0[-1]]
                HGVertForLookup[vert,0] = gridTotWatDepthVertex[vert,equalTo0[-1]]
                HWVertForLookup[vert,0] = wetTotWatDepthVertex[vert,equalTo0[-1]]
                cfVertForLookup[vert,0] = cfVertex[vert,equalTo0[-1]]
                cmfVertForLookup[vert,0] = cmfVertex[vert,equalTo0[-1]]
                
            else: # so if it never gets fully dry set everything to the value corresponding to the first surface elevations
            
                depthsVertForLookup[vert,0] = surfaceElevations[0]
                HGVertForLookup[vert,0] = gridTotWatDepthVertex[vert,0]
                HWVertForLookup[vert,0] = wetTotWatDepthVertex[vert,0]
                cfVertForLookup[vert,0] = cfVertex[vert,0]
                cmfVertForLookup[vert,0] = cmfVertex[vert,0]
                
            # now check for when phi == 1.0 and find exactly where that is
            
            equalTo1 = np.where(currPhiArray == 1.0)[0]
            
            if(len(equalTo1)!=0): # if 1.0 exists in the array
            
                depthsVertForLookup[vert,-1] = surfaceElevations[equalTo1[0]]
                HGVertForLookup[vert,-1] = gridTotWatDepthVertex[vert,equalTo1[0]]
                HWVertForLookup[vert,-1] = wetTotWatDepthVertex[vert,equalTo1[0]]
                cfVertForLookup[vert,-1] = cfVertex[vert,equalTo1[0]]
                cmfVertForLookup[vert,-1] = cmfVertex[vert,equalTo1[0]]
                
            else: # if there is nothing that is equal to 1 (so never gets fully wet, just set everything to correspind to the last surface elevation)
            
                depthsVertForLookup[vert,-1] = surfaceElevations[-1]
                HGVertForLookup[vert,-1] = gridTotWatDepthVertex[vert,-1]
                HWVertForLookup[vert,-1] = wetTotWatDepthVertex[vert,-1]
                cfVertForLookup[vert,-1] = cfVertex[vert,-1]
                cmfVertForLookup[vert,-1] = cmfVertex[vert,-1]
                
                
            # now for everything else
            
            for k in range(1,len(desiredPhiList)-1):
            
                desiredPhi = desiredPhiList[k]
                greaterThan = np.where(currPhiArray > desiredPhi)[0]
                
                if(len(greaterThan)==0): # so if nothing in the currPhiArray is greater than the desired phi
                
                # set everything to correspond to the last surface elevation
                    
                    depthsVertForLookup[vert,k] = surfaceElevations[-1]
                    HGVertForLookup[vert,k] = gridTotWatDepthVertex[vert,-1]
                    HWVertForLookup[vert,k] = wetTotWatDepthVertex[vert,-1]
                    cfVertForLookup[vert,k] = cfVertex[vert,-1]
                    cmfVertForLookup[vert,k] = cmfVertex[vert,-1]
                    
                elif(greaterThan[0] == 0): # so if the first currphi index is greater than the desired phi 
                
                # set everything to correspond to the first surfaceelevation
                
                    depthsVertForLookup[vert,k] = surfaceElevations[0]
                    HGVertForLookup[vert,k] = gridTotWatDepthVertex[vert,0]
                    HWVertForLookup[vert,k] = wetTotWatDepthVertex[vert,0]
                    cfVertForLookup[vert,k] = cfVertex[vert,0]
                    cmfVertForLookup[vert,k] = cmfVertex[vert,0]
        
                    
                else: # this is where we interpolate 
                    
                    greaterThan = greaterThan[0]
                    lessThan = greaterThan - 1
            
                    
                    depthsVertForLookup[vert,k] = (((desiredPhi - currPhiArray[lessThan])
                                                  /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                  *(surfaceElevations[greaterThan] - surfaceElevations[lessThan])
                                                  + (surfaceElevations[lessThan]))
                    
                    HGVertForLookup[vert,k] = (((desiredPhi - currPhiArray[lessThan])
                                                  /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                  *(gridTotWatDepthVertex[vert,greaterThan]
                                                    - gridTotWatDepthVertex[vert,lessThan])
                                                  + (gridTotWatDepthVertex[vert,lessThan]))
        
                    HWVertForLookup[vert,k] = (((desiredPhi - currPhiArray[lessThan])
                                                  /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                  *(wetTotWatDepthVertex[vert,greaterThan]
                                                    - wetTotWatDepthVertex[vert,lessThan])
                                                  + (wetTotWatDepthVertex[vert,lessThan]))
                    
                    cfVertForLookup[vert,k] = (((desiredPhi - currPhiArray[lessThan])
                                                  /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                  *(cfVertex[vert,greaterThan]
                                                    - cfVertex[vert,lessThan])
                                                  + (cfVertex[vert,lessThan]))
                    
                    cmfVertForLookup[vert,k] = (((desiredPhi - currPhiArray[lessThan])
                                                  /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                  *(cmfVertex[vert,greaterThan]
                                                    - cmfVertex[vert,lessThan])
                                                  + (cmfVertex[vert,lessThan]))
            
                

        # fill the vertices not contained in the subgrid region with -99999
        
        depthsVertForLookup[vertNotInSubgrid,:] = -99999
        HGVertForLookup[vertNotInSubgrid,:] = -99999
        HWVertForLookup[vertNotInSubgrid,:] = -99999
        cfVertForLookup[vertNotInSubgrid,:] = -99999
        cmfVertForLookup[vertNotInSubgrid,:] = -99999

        # deallocate arrays 
        wetFractionVertex = None
        gridTotWatDepthVertex = None
        wetTotWatDepthVertex = None
        cf = None
        cmf = None
        
        end = time.time()
        print('Reduction of partially wet vertices finished and took {} s'.format(end-start))
        
        #### CREATE THE NETCDF FILE
        ncFile = nc.Dataset(outputFilename, mode = 'w', format = 'NETCDF4')
        
        # create dimensions
        
        ncFile.createDimension('numEle',numEle) # element dimension
        ncFile.createDimension('numVert',3) # number of vertices per element
        ncFile.createDimension('numPhi',11) # number of possible phi values
        ncFile.createDimension('numSfcElevs',len(surfaceElevations)) # number of surface elevations
        ncFile.createDimension('numNode',numNode) # number of nodes in mesh
        ncFile.createDimension('oneDim',1) # just for a single constant to be saved
        
        # create variables
        # write variable dimensions transposed because FORTRAN will read them that way
        # only need to do it for the 3D arrays, handle 2 and 1D a different way.
        
        # wetAreaFraction set
        
        phiSet = ncFile.createVariable('phiSet',np.float32,
                                        'numPhi')
        
        # create array for depth in reduced vertex table
        wetFractionDepthVarVertex = ncFile.createVariable('wetFractionVertex',np.float32,
                                            ('numNode','numPhi'))
        
        # create array for depth in reduced element table
        wetFractionVarDepths = ncFile.createVariable('wetFractionDepths',np.float32,
                                                ('numPhi','numVert','numEle'))
        # elemental areas
        areaVar = ncFile.createVariable('area',np.float32,('numEle','numVert'))
        
        # elemental grid averaged total water depth
        totWatDepthVar = ncFile.createVariable('totWatDepth',np.float32,
                                                ('numPhi','numVert','numEle'))
        
        # # surface elevation array
        # surfaceElevationsVar = ncFile.createVariable('surfaceElevations',np.float32,'numSfcElevs')
        # minimum surface elevation
        minSurfElevVar = ncFile.createVariable('minSurfElev',np.float32,'oneDim')
        # maximum surface elevation
        maxSurfElevVar = ncFile.createVariable('maxSurfElev',np.float32,'oneDim')
        
        # vertex wet total water depth
        wetTotWatDepthVarVertex = ncFile.createVariable('wetTotWatDepthVertex',np.float32,
                                                          ('numNode','numPhi'))
        
        # vertex grid total water depth
        gridTotWatDepthVarVertex = ncFile.createVariable('gridTotWatDepthVertex',np.float32,
                                                          ('numNode','numPhi'))
        
        # vertex coefficient of friction level 0
        cfVarVertex = ncFile.createVariable('cfVertex',np.float32,
                                            ('numNode','numPhi'))
        
        # variables showing which elements and vertices are contained within
        # the subgrid area
                
        binaryElementListVariable = ncFile.createVariable('binaryElementList',np.int32,
                                                                  ('numEle'))
        binaryVertexListVariable = ncFile.createVariable('binaryVertexList',np.int32,
                                                                  ('numNode'))
        
        # do not use this in calculations anymore JLW 01302024
        # # write max Elevation
        
        # maxElevationEleVariable = ncFile.createVariable('maxElevationElement',np.float32,
        #                                               ('numEle'))
        
        # maxElevationVertexVariable = ncFile.createVariable('maxElevationVertex',np.float32,
        #                                               ('numNode'))
        # vertex coefficient of friction level 1
        cmfVarVertex = ncFile.createVariable('cmfVertex',np.float32,
                                              ('numNode','numPhi'))
        
        # elemental advection correction
        cadvVar = ncFile.createVariable('cadv',np.float32,
                                        ('numPhi','numVert','numEle'))
        
        
        phiSet[:] = desiredPhiList
        wetFractionDepthVarVertex[:,:] = depthsVertForLookup
        wetFractionVarDepths[:,:,:] = depthsEleForLookup
        areaVar[:,:] = area
        totWatDepthVar[:,:,:] = HEleForLookup
        # surfaceElevationsVar[:] = surfaceElevations
        minSurfElevVar[:] = minSurfElev
        maxSurfElevVar[:] = maxSurfElev
        wetTotWatDepthVarVertex[:,:] = HWVertForLookup
        gridTotWatDepthVarVertex[:,:] = HGVertForLookup
        cfVarVertex[:,:] = cfVertForLookup
        binaryElementListVariable[:] = binaryElementList
        binaryVertexListVariable[:] = binaryVertexList
        # do not use this in calculations anymore JLW 01302024
        # add max elevation cal
        # maxElevationEleVariable[:] = maxElevationEle
        # maxElevationVertexVariable[:] = maxElevationVertex
        cmfVarVertex[:,:] = cmfVertForLookup
        cadvVar[:,:,:] = cadvForLookup
        
        ncFile.close()
        

        ##################### SUBGRID CALULATOR NEW VERTEX SCHEME ###########################

    def calculateSubgridCorrection(controlFilename):

        import sys
        import numpy as np
        import time
        import matplotlib.pyplot as plt
        import re
        import netCDF4 as nc

        # first read in the control file
        
        controlFile = controlFilename
                
        with open(controlFile) as ctrF:
                    
            ctrF.readline()
            # change to shintaros r.strip with re to allow for spaces
            line = ctrF.readline().rstrip()
            line = re.split(' *= *',line)
            # get output file name
            outputFilename = line[1]
            line = ctrF.readline().rstrip()
            line = re.split(' *= *',line)
            # get mesh filename
            meshFilename = line[1]
            # read in mannings stuff
            line = ctrF.readline().rstrip()
            line = re.split(' *= *',line)
            defaultManning = line[1] # if true just use the manning table in the code
            if defaultManning == 'False': # otherwise we need to read a table in
                line = ctrF.readline().rstrip()
                line = re.split(' *= *',line)
                manningsnTableFilename = line[1] # get the mannings n table filename
            # now read in the elevation array
            line = ctrF.readline().rstrip()
            line = re.split(' *= *',line)
            minSurElev = float(line[1])
            line = ctrF.readline().rstrip()
            line = re.split(' *= *',line)
            maxSurElev = float(line[1])
            line = ctrF.readline().rstrip()
            line = re.split(' *= *',line)
            elevDisc = float(line[1])
            line = ctrF.readline().rstrip()
            line = re.split(' *= *',line)
            numDEMs = int(line[1])
            # get list of elevation datasets
            demFilenameList = []
            for i in range(numDEMs):
                line = ctrF.readline().rstrip()
                line = re.split(' *= *',line)
                demFilenameList.append(line[0])
            line = ctrF.readline().rstrip()
            line = re.split(' *= *',line)
            numLCs = int(line[1])
            # get list of landcover datasets
            landcoverFilenameList = []
            for i in range(numLCs):
                line = ctrF.readline().rstrip()
                line = re.split(' *= *',line)
                landcoverFilenameList.append(line[0])
        # get table to convert landcover values (possibly just make a function)
        if defaultManning:
                    
            landCoverToManning = {0:0.02, 2:0.15, 3:0.10, 4:0.05, 5:0.02,
                                6:0.037, 7:0.033, 8:0.034, 9:0.1, 10:0.11,
                                11:0.1, 12:0.05, 13:0.1, 14:0.048, 15:0.045,
                                16:0.1, 17:0.048, 18:0.045, 19:0.04,
                                20:0.09, 21:0.02, 22:0.015, 23:0.015, 
                                24:0.09, 25:0.01}
                        
        else:
                        
            landCoverToManning = subgridCalculatormain.readManning(manningsnTableFilename)

        landCoverValues = [0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25]

        # surface elevation array for calcuations
        surfaceElevations = np.round(np.arange(minSurElev,maxSurElev+elevDisc,
                                                elevDisc),2).astype('float32') 
        num_SfcElevs = len(surfaceElevations)

        # now read in the mesh 
        # meshFilename = r'C:\Users\jwoodruff\Desktop\subgrid_paper\simulations\SABv3-1000m\fort.14'
        mesh = subgridCalculatormain.readMesh(meshFilename)    
        meshConnectivity = mesh[1].triangles
        meshLon = np.asarray(mesh[0]['Longitude']).astype('float32')
        meshLat = np.asarray(mesh[0]['Latitude']).astype('float32')
        numNode = mesh[2]
        numEle = mesh[3]

        # for each vertex find what elements it is connected to

        # first mind the maximum number of elements connected to a vertex 
        # find most common vertex in connectivity table
        counts = np.bincount(meshConnectivity.flatten())
        commonVertex = np.argmax(counts)
        maxConnectedVertex = np.count_nonzero(meshConnectivity==commonVertex)
        # create array to hold all the vertex sub area data
        # 8 is used because each sub area will have 4 vertices so listing lon then lat
        vertexData = np.empty((numNode,maxConnectedVertex,8))
        vertexData[:,:,:] = np.nan

        # fill this vertex Data Array
        for i in range(numNode):
            # find connected elements (this is the slowest part)
            connectedElements = np.where(meshConnectivity==i)[0]
            # fill in vertex data
            for j in range(len(connectedElements)):
                # get vertices of subelement
                ele = connectedElements[j]
                # other vertices 
                otherVertices = meshConnectivity[ele,meshConnectivity[ele,:]!=i]
                # order vertices how you want
                # start with vertex in question
                nm0 = i
                nm1 = otherVertices[0]
                nm2 = otherVertices[1]
                # get lon and lat of vertices
                vertexNumbers = [nm0,nm1,nm2]
                vertexLon = meshLon[vertexNumbers]
                vertexLat = meshLat[vertexNumbers]
                # now get the centroid of the element
                centroidLon = (vertexLon[0]+vertexLon[1]+vertexLon[2])/3
                centroidLat = (vertexLat[0]+vertexLat[1]+vertexLat[2])/3
                # get mid point of each vertex connected to vertex of interest
                midPointLon1 = (meshLon[nm0]+meshLon[nm1])/2
                midPointLon2 = (meshLon[nm0]+meshLon[nm2])/2
                midPointLat1 = (meshLat[nm0]+meshLat[nm1])/2
                midPointLat2 = (meshLat[nm0]+meshLat[nm2])/2
                # now add this data to vertex array
                subAreaPerimeter = np.array((meshLon[nm0],midPointLon1,centroidLon,midPointLon2,
                                            meshLat[nm0],midPointLat1,centroidLat,midPointLat2))
                vertexData[i,j,:] = subAreaPerimeter

        # get an array of max and min vertex area coordinates
        vertexAreaMinLon = np.nanmin(vertexData[:,:,:4],axis=(1,2))
        vertexAreaMaxLon = np.nanmax(vertexData[:,:,:4],axis=(1,2))
        vertexAreaMinLat = np.nanmin(vertexData[:,:,4:],axis=(1,2))
        vertexAreaMaxLat = np.nanmax(vertexData[:,:,4:],axis=(1,2))

        # preallocate arrays to store vertex subgrid data
        wetFractionVertex = np.empty((numNode,num_SfcElevs))
        wetFractionVertex[:,:] = np.nan
        totWatDepthVertex = np.empty(wetFractionVertex.shape)
        totWatDepthVertex[:,:] = np.nan
        wetTotWatDepthVertex = np.empty(wetFractionVertex.shape)
        wetTotWatDepthVertex[:,:] = np.nan
        cfVertex = np.empty(wetFractionVertex.shape)
        cfVertex[:,:] = np.nan
        cmfVertex = np.empty(wetFractionVertex.shape)
        cmfVertex[:,:] = np.nan
        cadvVertex = np.empty(wetFractionVertex.shape)
        cadvVertex[:,:] = np.nan
        # create a list of vertices contained within the subgrid area
        vertexUseList = np.ones(numNode,dtype=bool)
        # keep track of total calc time
        startTotal = time.time()
        for i in range(len(demFilenameList)):
            # reading in DEM again
            # all variables the same as before
            elevationData = subgridCalculatormain.importDEMv2(demFilenameList[i])
            landcoverData = subgridCalculatormain.importDEMv2(landcoverFilenameList[i])
            # get data out of dems
            bathyTopo = elevationData[0].astype('float32')
            lonRes = elevationData[1]
            latRes = -1*elevationData[2]
            lon = elevationData[3]
            lat = elevationData[4]
            elevationData = None # deallocate 
            manningsn = landcoverData[0].astype('float32') # array of mannings n values
            landcoverData = None # deallocate
            
            for value in landCoverValues:
                    
                manningsn[manningsn == value] = landCoverToManning[value]
                    
            # set nan values to 0.02
                    
            manningsn[np.isnan(manningsn)] = 0.02
            # create a buffer for containment
            # bufferCells = 3
            bufferLon = lonRes# *bufferCells
            bufferLat = latRes# *bufferCells
            # get bounds of raster data
            demMinLon = np.min(lon)
            demMaxLon = np.max(lon)
            demMinLat = np.min(lat)
            demMaxLat = np.max(lat)

            # see what vertices lie within this dem
            minLonWithin = np.array(vertexAreaMinLon>demMinLon)
            maxLonWithin = np.array(vertexAreaMaxLon<demMaxLon)
            minLatWithin = np.array(vertexAreaMinLat>demMinLat)
            maxLatWithin = np.array(vertexAreaMaxLat<demMaxLat)
            # get all within and only use vertices that have not already been used
            allWithin = minLonWithin*maxLonWithin*minLatWithin*maxLatWithin*vertexUseList
            idxAllWithin = np.where(allWithin)[0]
            # update vertex use list 
            vertexUseList[idxAllWithin] = False

            # now loop through contained vertices and perform calculations
            # for vertex in idxAllWithin:
            for j in range(len(idxAllWithin)):
                start = time.time()
                # find how many connected elements
                conElementCount = np.count_nonzero(~np.isnan(vertexData[idxAllWithin[j],:,0])) 
                # create array to hold vertex area data 
                vertexSubArea = np.zeros((conElementCount,1))
                # temporarily allocate for each subarea variable
                tempwetFractionData = np.zeros((conElementCount,num_SfcElevs))
                # just set the rest equal to this because we are just preallocating 0 arrays
                temptotWatDepthData = np.zeros((conElementCount,num_SfcElevs)) 
                tempwetTotWatDepthData = np.zeros((conElementCount,num_SfcElevs)) 
                tempcfData = np.zeros((conElementCount,num_SfcElevs)) 
                tempcmfData = np.zeros((conElementCount,num_SfcElevs)) 
                tempcadvData = np.zeros((conElementCount,num_SfcElevs)) 
                # now loop through connected elements
                for k in range(conElementCount):
                    # get the sub area perimeter for the particular element
                    subAreaPerimeterLon = vertexData[idxAllWithin[j],k,:4]
                    subAreaPerimeterLat = vertexData[idxAllWithin[j],k,4:]
                    # cut down dem and landcover to each sub area
                    # get locations of bounds
                    minLonDEMWithin = lon > np.min(subAreaPerimeterLon)
                    maxLonDEMWithin = lon < np.max(subAreaPerimeterLon)
                    lonWithinIdx = np.where(minLonDEMWithin*maxLonDEMWithin)[0]
                    minCol = np.min(lonWithinIdx)#  - bufferCells # create a cell buffer
                    maxCol = np.max(lonWithinIdx)#  + bufferCells
                    demLonCut = lon[minCol:maxCol]
                    minLatDEMWithin = lat > np.min(subAreaPerimeterLat)
                    maxLatDEMWithin = lat < np.max(subAreaPerimeterLat)
                    latWithinIdx = np.where(minLatDEMWithin*maxLatDEMWithin)[0]
                    minRow = np.min(latWithinIdx)#  - bufferCells # create a 3 cell buffer
                    maxRow = np.max(latWithinIdx)#  + bufferCells
                    demLatCut = lat[minRow:maxRow]
                    demBathyTopoCut = bathyTopo[minRow:maxRow,minCol:maxCol]
                    manningsnCut = manningsn[minRow:maxRow,minCol:maxCol]
                    lonGrid,latGrid = np.meshgrid(demLonCut,demLatCut)

                    # split into 2 triangles
                        
                    triLon0 = subAreaPerimeterLon[:3]
                    triLat0 = subAreaPerimeterLat[:3]
                            
                    # convert to meters
                    tri0Meters = subgridCalculatormain.projectMeshToMercator(triLat0,
                                                                            triLon0)
                    tri0Area = subgridCalculatormain.triarea(tri0Meters[0][0], 
                                                            tri0Meters[1][0],
                                                            tri0Meters[0][1], 
                                                            tri0Meters[1][1],
                                                            tri0Meters[0][2], 
                                                            tri0Meters[1][2])

                            
                    insideTri0 = subgridCalculatormain.isInside(triLon0[0], triLat0[0],
                                                                triLon0[1], triLat0[1],
                                                                triLon0[2], triLat0[2],
                                                                lonGrid, latGrid, 0.00000001)
                            
                    triLon1 = subAreaPerimeterLon[[0,2,3]]
                    triLat1 = subAreaPerimeterLat[[0,2,3]]
                    tri1Meters = subgridCalculatormain.projectMeshToMercator(triLat1,
                                                                            triLon1)

                            
                    tri1Area = subgridCalculatormain.triarea(tri1Meters[0][0], 
                                                            tri1Meters[1][0],
                                                            tri1Meters[0][1], 
                                                            tri1Meters[1][1],
                                                            tri1Meters[0][2], 
                                                            tri1Meters[1][2])
                            
                            
                    insideTri1 = subgridCalculatormain.isInside(triLon1[0], triLat1[0],
                                                                triLon1[1], triLat1[1],
                                                                triLon1[2], triLat1[2],
                                                                lonGrid, latGrid, 0.00000001)
                            
                    # now combine the two triangles and find the points inside
                            
                    insideSubElement = np.logical_or(insideTri0,insideTri1)
                            
                    # count the number of subgrid cells within the subelement
                            
                    cellsInSubElement = np.count_nonzero(insideSubElement)

                    # if there are no cells within the element the DEM is too coarse
                    # you must decrease the DEM resolution in this area
                    if cellsInSubElement == 0:

                        sys.exit('DEM {0} resolution too coarse!'.format(i))
                            
                    # get just he bathy topo inside the sub element 
                            
                    bathyTopoInsideSubElement = demBathyTopoCut*insideSubElement
                            
                    # set 0 values to nan for calculations
                            
                    bathyTopoInsideSubElement[bathyTopoInsideSubElement==0] = np.nan
                            
                    # get area of sub element
                    vertexSubArea[k] = tri0Area + tri1Area # used for area weighting later
                    # remove cells not inside sub element which will flatten the array
                    bathyTopoInsideSubElementNoNaN = bathyTopoInsideSubElement[~np.isnan(bathyTopoInsideSubElement)]
                    manningsnCutNoNaN = manningsnCut[~np.isnan(bathyTopoInsideSubElement)]

                    # get the total water depth at each surface elevation
                            
                    temptotWatDepth =  surfaceElevations[:,None] - bathyTopoInsideSubElementNoNaN
                            
                    # count the number of wet cells
                            
                    wetCellsInSubArea = temptotWatDepth > 0.0001
                            
                    wetCellsInSubAreaCount = np.sum(wetCellsInSubArea,axis=1).astype('float32')

                    # now set tot water depth of dry cells to nan
                            
                    temptotWatDepth[temptotWatDepth < 0.0001] = np.nan
            
                    # add to wet frac array
                            
                    tempwetFractionData[k,:] = wetCellsInSubAreaCount/cellsInSubElement

                    # add to total water depth array                       
                    temptotWatDepthData[k,:] = np.nansum(temptotWatDepth,axis=1)/cellsInSubElement

                    # get wet total water depth and coefficient of friction
                            
                    # find the mannings for only wet areas then 0 the rest for 
                    # use in calculations 
                            
                    manningsnCutNoNaNWet = manningsnCutNoNaN * wetCellsInSubArea
                            
                    tempcf = (9.81*manningsnCutNoNaNWet**2)/(temptotWatDepth**(1/3))
                    # set 0 tempcf to nan to prevent 0 divide
                    tempcf[tempcf==0] = np.nan

                    # make wet cells in sub area count nan when == 0 so we don't get divide by 0 issue
                    wetCellsInSubAreaCount[wetCellsInSubAreaCount==0.0] = np.nan
                    # calculate wet total water depth
                    tempwetTotWatDepthData[k,:] = np.nansum(temptotWatDepth,axis=1)/wetCellsInSubAreaCount
                    # calculate bottom friction coefficient
                    tempcfData[k,:]= np.nansum(tempcf,axis=1)/wetCellsInSubAreaCount
                    # calculate rv for use in advection and bottom friction correction
                    rv = tempwetTotWatDepthData[k,:]/(np.nansum((temptotWatDepth**(3/2))*(tempcf**(-1/2)),axis=1)/wetCellsInSubAreaCount)
                    # calcualte advection correction
                    tempcadvData[k,:] = (1/tempwetTotWatDepthData[k,:])*(np.nansum(temptotWatDepth**2/tempcf,axis=1)/wetCellsInSubAreaCount)*rv**2
                    # calculate bottom friction correction
                    tempcmfData[k,:] = tempwetTotWatDepthData[k,:]*rv**2 # this is correct I need <H>W * Rv**2
                    # now fill the nans from this calculation which represent depths were nothing was wet
                    # set wet total water depth to 0
                    tempwetTotWatDepthData[k,np.isnan(tempwetTotWatDepthData[k,:])] = 0.0
                    # set nan values to cf calculated from mean mannings n and 8 cm of water
                    tempcfData[k,np.isnan(tempcfData[k,:])] = 9.81*np.mean(manningsnCutNoNaN)**2/(0.08**(1/3))
                    tempcmfData[k,np.isnan(tempcmfData[k,:])] = 9.81*np.mean(manningsnCutNoNaN)**2/(0.08**(1/3))
                    # set advection correction equal to 1.0
                    tempcadvData[k,np.isnan(tempcadvData[k,:])] = 1.0
                
                # once the sub elements surrounding each vertex have been looped through, put all of he data on the elements
                areaTotalVertex = np.sum(vertexSubArea[:,0])
        
                # check = np.sum(tempwetFractionData*vertexSubArea,axis=0)/areaTotalVertex
                wetFractionVertex[idxAllWithin[j],:] = np.sum(tempwetFractionData*vertexSubArea,axis=0)/areaTotalVertex
                totWatDepthVertex[idxAllWithin[j],:] = np.sum(temptotWatDepthData*vertexSubArea,axis=0)/areaTotalVertex
                wetTotWatDepthVertex[idxAllWithin[j],:] = np.sum(tempwetTotWatDepthData*vertexSubArea,axis=0)/areaTotalVertex
                cfVertex[idxAllWithin[j],:] = np.sum(tempcfData*vertexSubArea,axis=0)/areaTotalVertex
                cmfVertex[idxAllWithin[j],:] = np.sum(tempcmfData*vertexSubArea,axis=0)/areaTotalVertex
                cadvVertex[idxAllWithin[j],:] = np.sum(tempcadvData*vertexSubArea,axis=0)/areaTotalVertex
                # finish vertex loop and print time
                end = time.time()
                print('Finish vertex {} in DEM {} took : {} s'.format(idxAllWithin[j],i,end-start))

        # set minimums for bottom friction
        cfVertex[cfVertex<0.0025] = 0.0025
        cmfVertex[cmfVertex<0.0025] = 0.0025
        # total time
        endTotal = time.time()
        print('All calulations took {} s'.format(endTotal-startTotal))


        # now I need to condense the lookup tables to only have 11 values corresponding to phi=0 to phi=1
        start = time.time() 
        desiredPhiList = np.linspace(0,1,11)
        depthsVertForLookup = np.ones((len(wetFractionVertex[:]),11))*-99999
        HGVertForLookup = np.ones((len(wetFractionVertex[:]),11))*-99999
        HWVertForLookup = np.ones((len(wetFractionVertex[:]),11))*-99999
        cfVertForLookup = np.ones((len(wetFractionVertex[:]),11))*-99999
        cmfVertForLookup = np.ones((len(wetFractionVertex[:]),11))*-99999
        cadvVertForLookup = np.ones((len(wetFractionVertex[:]),11))*1.0

        # only loop through Nodes in subgrid area 
        # get list of which vertices are inside the subgrid area
        vertsInSubArea = np.where(vertexUseList==False)[0] # confusing that I use false but it was for multiplication earlier

        for vert in vertsInSubArea:

            currPhiArray = wetFractionVertex[vert,:]
                    
            # make sure that the phi array also gets fully wet and then proceed
            # otherwise just skip
                    
            # for phi == 0 you want to find exactly where that is in the currPhiArray
            equalTo0 = np.where(currPhiArray == 0.0)[0]
                    
            if(len(equalTo0)!=0): # if 0.0 exists in the array
                    
                depthsVertForLookup[vert,0] = surfaceElevations[equalTo0[-1]]
                HGVertForLookup[vert,0] = totWatDepthVertex[vert,equalTo0[-1]]
                HWVertForLookup[vert,0] = wetTotWatDepthVertex[vert,equalTo0[-1]]
                cfVertForLookup[vert,0] = cfVertex[vert,equalTo0[-1]]
                cmfVertForLookup[vert,0] = cmfVertex[vert,equalTo0[-1]]
                cadvVertForLookup[vert,0] = cadvVertex[vert,equalTo0[-1]]
                        
            else: # so if it never gets fully dry set everything to the value corresponding to the first surface elevations
                    
                depthsVertForLookup[vert,0] = surfaceElevations[0]
                HGVertForLookup[vert,0] = totWatDepthVertex[vert,0]
                HWVertForLookup[vert,0] = wetTotWatDepthVertex[vert,0]
                cfVertForLookup[vert,0] = cfVertex[vert,0]
                cmfVertForLookup[vert,0] = cmfVertex[vert,0]
                cadvVertForLookup[vert,0] = cadvVertex[vert,0]
                        
            # now check for when phi == 1.0 and find exactly where that is
                    
            equalTo1 = np.where(currPhiArray == 1.0)[0]
                    
            if(len(equalTo1)!=0): # if 1.0 exists in the array
                    
                depthsVertForLookup[vert,-1] = surfaceElevations[equalTo1[0]]
                HGVertForLookup[vert,-1] = totWatDepthVertex[vert,equalTo1[0]]
                HWVertForLookup[vert,-1] = wetTotWatDepthVertex[vert,equalTo1[0]]
                cfVertForLookup[vert,-1] = cfVertex[vert,equalTo1[0]]
                cmfVertForLookup[vert,-1] = cmfVertex[vert,equalTo1[0]]
                cadvVertForLookup[vert,-1] = cadvVertex[vert,equalTo1[0]]
                        
            else: # if there is nothing that is equal to 1 (so never gets fully wet, just set everything to correspind to the last surface elevation)
                    
                depthsVertForLookup[vert,-1] = surfaceElevations[-1]
                HGVertForLookup[vert,-1] = totWatDepthVertex[vert,-1]
                HWVertForLookup[vert,-1] = wetTotWatDepthVertex[vert,-1]
                cfVertForLookup[vert,-1] = cfVertex[vert,-1]
                cmfVertForLookup[vert,-1] = cmfVertex[vert,-1]
                cadvVertForLookup[vert,-1] = cadvVertex[vert,-1]
                        
                        
            # now for everything else
                    
            for k in range(1,len(desiredPhiList)-1):
                    
                desiredPhi = desiredPhiList[k]
                greaterThan = np.where(currPhiArray > desiredPhi)[0]
                        
                if(len(greaterThan)==0): # so if nothing in the currPhiArray is greater than the desired phi
                        
                # set everything to correspond to the last surface elevation
                            
                    depthsVertForLookup[vert,k] = surfaceElevations[-1]
                    HGVertForLookup[vert,k] = totWatDepthVertex[vert,-1]
                    HWVertForLookup[vert,k] = wetTotWatDepthVertex[vert,-1]
                    cfVertForLookup[vert,k] = cfVertex[vert,-1]
                    cmfVertForLookup[vert,k] = cmfVertex[vert,-1]
                    cadvVertForLookup[vert,k] = cadvVertex[vert,-1]
                            
                elif(greaterThan[0] == 0): # so if the first currphi index is greater than the desired phi 
                        
                # set everything to correspond to the first surfaceelevation
                        
                    depthsVertForLookup[vert,k] = surfaceElevations[0]
                    HGVertForLookup[vert,k] = totWatDepthVertex[vert,0]
                    HWVertForLookup[vert,k] = wetTotWatDepthVertex[vert,0]
                    cfVertForLookup[vert,k] = cfVertex[vert,0]
                    cmfVertForLookup[vert,k] = cmfVertex[vert,0]
                    cadvVertForLookup[vert,k] = cadvVertex[vert,0]
                
                            
                else: # this is where we interpolate 
                            
                    greaterThan = greaterThan[0]
                    lessThan = greaterThan - 1
            
                    depthsVertForLookup[vert,k] = (((desiredPhi - currPhiArray[lessThan])
                                                /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                *(surfaceElevations[greaterThan] - surfaceElevations[lessThan])
                                                + (surfaceElevations[lessThan]))
                    HGVertForLookup[vert,k] = (((desiredPhi - currPhiArray[lessThan])
                                                /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                *(totWatDepthVertex[vert,greaterThan]
                                                    - totWatDepthVertex[vert,lessThan])
                                                + (totWatDepthVertex[vert,lessThan]))
                    HWVertForLookup[vert,k] = (((desiredPhi - currPhiArray[lessThan])
                                                /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                *(wetTotWatDepthVertex[vert,greaterThan]
                                                - wetTotWatDepthVertex[vert,lessThan])
                                                + (wetTotWatDepthVertex[vert,lessThan]))
                    cfVertForLookup[vert,k] = (((desiredPhi - currPhiArray[lessThan])
                                                /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                *(cfVertex[vert,greaterThan]
                                                    - cfVertex[vert,lessThan])
                                                + (cfVertex[vert,lessThan]))
                    cmfVertForLookup[vert,k] = (((desiredPhi - currPhiArray[lessThan])
                                                /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                *(cmfVertex[vert,greaterThan]
                                                - cmfVertex[vert,lessThan])
                                                + (cmfVertex[vert,lessThan]))
                    cadvVertForLookup[vert,k] = (((desiredPhi - currPhiArray[lessThan])
                                                /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                *(cadvVertex[vert,greaterThan]
                                                - cadvVertex[vert,lessThan])
                                                + (cadvVertex[vert,lessThan]))
                
        end = time.time()
        print('Reduction of partially wet vertices finished and took {} s'.format(end-start))

        #### CREATE THE NETCDF FILE
        ncFile = nc.Dataset(outputFilename, mode = 'w', format = 'NETCDF4')
                
        # create dimensions
        ncFile.createDimension('numPhi',11) # number of possible phi values
        ncFile.createDimension('numNode',numNode) # number of nodes in mesh
                
        # create variables
        # write variable dimensions transposed because FORTRAN will read them that way
        # only need to do it for the 3D arrays, handle 2 and 1D a different way.
                
        # wetAreaFraction set
                
        phiSet = ncFile.createVariable('phiSet',np.float32,
                                                'numPhi')
                
        # create array for depth in reduced vertex table
        wetFractionDepthVarVertex = ncFile.createVariable('wetFractionVertex',np.float32,
                                                    ('numNode','numPhi'))

        # vertex wet total water depth
        wetTotWatDepthVarVertex = ncFile.createVariable('wetTotWatDepthVertex',np.float32,
                                                                ('numNode','numPhi'))
                
        # vertex grid total water depth
        gridTotWatDepthVarVertex = ncFile.createVariable('gridTotWatDepthVertex',np.float32,
                                                                ('numNode','numPhi'))
                
        # vertex coefficient of friction level 0
        cfVarVertex = ncFile.createVariable('cfVertex',np.float32,
                                                    ('numNode','numPhi'))
                
        # variables showing which vertices are contained within
        # the subgrid area
        binaryVertexListVariable = ncFile.createVariable('binaryVertexList',np.int32,
                                                                        ('numNode'))
                
        # vertex coefficient of friction level 1
        cmfVarVertex = ncFile.createVariable('cmfVertex',np.float32,
                                                    ('numNode','numPhi'))
        # vertex advection correction
        cadvVarVertex = ncFile.createVariable('cadvVertex',np.float32,
                                                    ('numNode','numPhi'))
        # create a binary list of vertices within the subgrid
        # if inside set to 1 otherwise 0
        binaryVertexList = np.zeros(len(vertexUseList)).astype('int')
        binaryVertexList[vertexUseList==False] = 1

        phiSet[:] = desiredPhiList
        wetFractionDepthVarVertex[:,:] = depthsVertForLookup
        wetTotWatDepthVarVertex[:,:] = HWVertForLookup
        gridTotWatDepthVarVertex[:,:] = HGVertForLookup
        cfVarVertex[:,:] = cfVertForLookup
        binaryVertexListVariable[:] = binaryVertexList
        cmfVarVertex[:,:] = cmfVertForLookup
        cadvVarVertex[:,:] = cadvVertForLookup
                
        ncFile.close()