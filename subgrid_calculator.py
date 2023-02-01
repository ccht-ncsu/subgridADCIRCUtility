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

##################### SUBGRID CPU CALULATOR V3 ################################

    def calculateSubgridCorrectionv3(controlFilename):
        
        import sys
        import netCDF4 as nc
        import numpy as np
        import time
        import matplotlib.pyplot as plt
        # use re to read control file in per Shintaro recommendation
        import re
        # import subgrid_calculator as sc
        
        # I want to really dig into speed on this code and figure out how to make the
        # the code fast
        
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
        
        # now create a surface elevation array
        
        ds = 0.2 # will want to experiment with this and speed
        # surface elevation array for caoluations
        surfaceElevations = np.round(np.arange(-5,5+ds,ds),2).astype('float32') 
        
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
        
        maxElevationEle = np.zeros(numEle).astype(np.float32) # find highest elevation in each element for use in variable phi
        
        maxElevationSubEle = np.zeros((numEle,3)).astype(np.float32) # find the highest elevation in each subElement
        
        # now fill the rows and columns of non subgrid vertices and elements with -99999
        
        wetFraction[np.where(binaryElementList == 0),:,:] = -99999
        area[np.where(binaryElementList == 0)] = -99999
        totWatDepth[np.where(binaryElementList == 0),:,:] = -99999
        wetTotWatDepth[np.where(binaryElementList == 0),:,:] = -99999
        cf[np.where(binaryElementList == 0),:,:] = -99999
        
        # fill max elevation
        maxElevationEle[np.where(binaryElementList == 0)] = -99999
        maxElevationSubEle[np.where(binaryElementList == 0)] = -99999
        
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
                    
                    # get maximum elevation inside the sub element
                    
                    maxElevationSubEle[ele,j] = np.nanmax(bathyTopoInsideSubElement)
                    
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
                    
                    wetCellsInSubArea = temptotWatDepth > 0.001
                    
                    wetCellsInSubAreaCount = np.sum(wetCellsInSubArea,axis=1)
                    
                    # now set tot water depth of dry cells to nan
                    
                    temptotWatDepth[temptotWatDepth < 0.001] = np.nan
                    
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
                    
                    cmf[ele,j,:] = wetTotWatDepth[ele,j,:]*rv[ele,j,:]**2 # this is corrrect I need <H>W * Rv**2
                    # cmf[ele,j,:] = totWatDepth[ele,j,:]*rv[ele,j,:]**2 # this is incorrect
                    
                
                # get the maximum elevation inside the element
                maxElevationEle[ele] = np.max(maxElevationSubEle[ele,:])
            
                                        
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
        
        # now put sub-elemental quantities on the vertices
        
        wetFractionVertex = np.zeros((numNode,num_SfcElevs))
        wetTotWatDepthVertex = np.zeros((numNode,num_SfcElevs))
        gridTotWatDepthVertex = np.zeros((numNode,num_SfcElevs))
        vertexArea = np.zeros((numNode,num_SfcElevs))
        cfVertex = np.zeros((numNode,num_SfcElevs))
        cmfVertex = np.zeros((numNode,num_SfcElevs))
        # cadvVertex = np.zeros((numNode,num_SfcElevs))
        maxElevationVertex = np.ones(numNode)*-99999
        
        # now fill non subgrid spaces with -99999
        
        wetFractionVertex[np.where(binaryVertexList == 0),:] = -99999
        wetTotWatDepthVertex[np.where(binaryVertexList == 0),:] = -99999
        gridTotWatDepthVertex[np.where(binaryVertexList == 0),:] = -99999
        vertexArea[np.where(binaryVertexList == 0),:] = -99999
        cfVertex[np.where(binaryVertexList == 0),:] = -99999
        maxElevationVertex[np.where(binaryVertexList == 0)] = -99999
        
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
                
                # add max elevation
                if(maxElevationVertex[nm0] < maxElevationSubEle[i,0]):
                    
                    maxElevationVertex[nm0] = maxElevationSubEle[i,0]
            
                if(maxElevationVertex[nm1] < maxElevationSubEle[i,1]):
                    
                    maxElevationVertex[nm1] = maxElevationSubEle[i,1]
                    
                if(maxElevationVertex[nm2] < maxElevationSubEle[i,2]):
                    
                    maxElevationVertex[nm2] = maxElevationSubEle[i,2]
                
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
                
                # cadv0 = cadv[i,0,:]
                # cadv1 = cadv[i,1,:]
                # cadv2 = cadv[i,2,:]
            
                cmfVertex[nm0,:] += cmf0 * area[i,0]
                cmfVertex[nm1,:] += cmf1 * area[i,1]
                cmfVertex[nm2,:] += cmf2 * area[i,2]
            
                # cadvVertex[nm0,:] += cadv0 * area[i,0]
                # cadvVertex[nm1,:] += cadv1 * area[i,1]
                # cadvVertex[nm2,:] += cadv2 * area[i,2]
        
        # now average all of these by the vertex areas
        wetFractionVertex[np.where(binaryVertexList == 1)] = wetFractionVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        wetTotWatDepthVertex[np.where(binaryVertexList == 1)] = wetTotWatDepthVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        gridTotWatDepthVertex[np.where(binaryVertexList == 1)] = gridTotWatDepthVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        cfVertex[np.where(binaryVertexList == 1)] = cfVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        cmfVertex[np.where(binaryVertexList == 1)] = cmfVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        # cadvVertex[np.where(binaryVertexList == 1)] = cadvVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        
        # now we need to check to see if there are any vertices in the subgrid
        # but not connected to any elements in the subgrid 
        
        for i in range(numNode):
            
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
                    cmfVertex[i,:] = -99999
                    # cadvVertex[i,:] = -99999
                        
            # otherwise if vertex is outside subgrid area
            # just make everything -99999
            else:
                
                # change the vertex variables to be -99999
                wetFractionVertex[i,:] = -99999
                wetTotWatDepthVertex[i,:] = -99999
                gridTotWatDepthVertex[i,:] = -99999
                cfVertex[i,:] = -99999
                cmfVertex[i,:] = -99999
                # cadvVertex[i,:] = -99999
                
        end = time.time()
        print('Averaged to vertices took {} s'.format(end-start))
        # now we need to simplify this lookup table
        start = time.time()
        desiredPhiList = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        
        # transpose the dimensions of the elemental subgrid arrays for writing 
        # netCDF
        
        wetFraction = wetFraction.T
        totWatDepth = totWatDepth.T
        cadv = cadv.T
        
        
        # create empty arrays for the reduced tables
        
        depthsEleForLookup = np.zeros((11,3,len(wetFraction[0,0,:])))
        HEleForLookup = np.zeros((11,3,len(wetFraction[0,0,:])))
        cadvForLookup = np.zeros((11,3,len(wetFraction[0,0,:])))
                            
        
        # now loop through the elements 
        
        for i in range(len(wetFraction[0,0,:])):
            
            for j in range(3):
                
                currPhiArray = wetFraction[:,j,i]
                
                # check to see if the element is always fully wet
                
                if(all(currPhiArray == 1.0)):
                    
                    # set all to really low depths so that in the code you know to set at 1
                
                    depthsEleForLookup[:,j,i] = np.min(surfaceElevations)
                    
                    # set the total water depth to the first value so that we have a reference
                    # point to add zeta to
                    
                    HEleForLookup[:,j,i] = totWatDepth[0,j,i]
                    
                    # set cadvection to the last value in the lookup table
                    
                    cadvForLookup[:,j,i] = cadv[-1,j,i]
                    
                # check to see if the element is always fully dry
                    
                elif(all(currPhiArray == 0.0)):
                    
                    # set all to really high dpeths so that in code you know to set to 0.0
                    
                    depthsEleForLookup[:,j,i] = np.max(surfaceElevations) 
                    HEleForLookup[:,j,i] = 0.0
                    
                    # set cadv to the first value of the lookup table
                    
                    cadvForLookup[:,j,i] = cadv[0,j,i]
        
                # if the element is partially wet for all values
                            
                elif(all(currPhiArray > 0.0) and all(currPhiArray < 1.0)):
                    
                    # set to always dry
                    
                    depthsEleForLookup[:,j,i] = np.max(surfaceElevations) 
                    HEleForLookup[:,j,i] = 0.0
                    cadvForLookup[:,j,i] = cadv[0,j,i]
                            
                # if partially wet to fully wet
                
                elif(all(currPhiArray > 0.0) and any(currPhiArray == 1.0)):
                    
                    # set to always wet
                    
                    depthsEleForLookup[:,j,i] = np.min(surfaceElevations)
                    
                    # set the total water depth to the first value so that we have a reference
                    # point to add zeta to
                    
                    HEleForLookup[:,j,i] = totWatDepth[0,j,i]
                    
                    # set cadvection to the last value in the lookup table
                    
                    cadvForLookup[:,j,i] = cadv[-1,j,i]
                            
                # if fully dry to partially wet
                
                elif(all(currPhiArray < 1.0) and any(currPhiArray == 0.0)):
                    
                    # set to always dry
                    
                    depthsEleForLookup[:,j,i] = np.max(surfaceElevations)  
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
                            depthsEleForLookup[k,j,i] = surfaceElevations[equalTo]
                            HEleForLookup[k,j,i] = totWatDepth[equalTo,j,i]
                            cadvForLookup[k,j,i] = cadv[equalTo,j,i]
                            
                        # for phi == 1.0 you want exactly where that is in the currPhiArray
                            
                        elif(desiredPhi == 1.0):
                            
                            equalTo = np.where(currPhiArray == desiredPhi)[0][0]
                            depthsEleForLookup[k,j,i] = surfaceElevations[equalTo]
                            HEleForLookup[k,j,i] = totWatDepth[equalTo,j,i]
                            cadvForLookup[k,j,i] = cadv[equalTo,j,i]
                            
                        # now look for the points in between 0 and 1
                            
                        else:
                            
                            greaterThan = np.where(currPhiArray > desiredPhi)[0][0]
                            lessThan = np.where(currPhiArray < desiredPhi)[0][-1]
                            
                            depthsEleForLookup[k,j,i] = (((desiredPhi - currPhiArray[lessThan])
                                                          /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                          *(surfaceElevations[greaterThan] - surfaceElevations[lessThan])
                                                          + (surfaceElevations[lessThan]))
                            
                            HEleForLookup[k,j,i] = (((desiredPhi - currPhiArray[lessThan])
                                                          /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                          *(totWatDepth[greaterThan,j,i] - totWatDepth[lessThan,j,i])
                                                          + (totWatDepth[lessThan,j,i]))
                            
                            cadvForLookup[k,j,i] = (((desiredPhi - currPhiArray[lessThan])
                                                          /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                          *(cadv[greaterThan,j,i] - cadv[lessThan,j,i])
                                                          + (cadv[lessThan,j,i]))
                            
                            
                                     
            # print('Element {0} completed'.format(i))
            
        # fill the elements that are not contained in the subgrid region with -99999
        depthsEleForLookup[:,:,np.where(binaryElementList==0)[0]] = -99999
        HEleForLookup[:,:,np.where(binaryElementList==0)[0]] = -99999
        cadvForLookup[:,:,np.where(binaryElementList==0)[0]] = -99999
        end = time.time()
        print('Reducing Lookuptable took {} s'.format(end-start))
        ncFile = nc.Dataset(outputFilename, mode = 'w', format = 'NETCDF4')
        
        # create dimensions
        
        ncFile.createDimension('numEle',numEle) # element dimension
        ncFile.createDimension('numVert',3) # number of vertices per element
        ncFile.createDimension('numPhi',11) # number of possible phi values
        ncFile.createDimension('numSfcElevs',len(surfaceElevations)) # number of surface elevations
        ncFile.createDimension('numNode',numNode) # number of nodes in mesh
        
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
        # vertex coefficient of friction level 1
        cmfVarVertex = ncFile.createVariable('cmfVertex',np.float32,
                                              ('numNode','numSfcElevs'))
        
        # elemental advection correction
        cadvVar = ncFile.createVariable('cadv',np.float32,
                                        ('numPhi','numVert','numEle'))
        
        
        phiSet[:] = desiredPhiList
        wetFractionVarVertex[:,:] = wetFractionVertex
        wetFractionVarDepths[:,:,:] = depthsEleForLookup
        areaVar[:,:] = area
        totWatDepthVar[:,:,:] = HEleForLookup
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
        cmfVarVertex[:,:] = cmfVertex
        cadvVar[:,:,:] = cadvForLookup
        
        ncFile.close()
        
##################### SUBGRID GPU CALULATOR V3 ################################
# Now use the v3 code and add GPU functionality

    def calculateSubgridCorrectionv3GPU(controlFilename):
        
        import sys
        import netCDF4 as nc
        import numpy as np
        import time
        import matplotlib.pyplot as plt
        import cupy as cp
        # use re to read control file in per Shintaro recommendation
        import re
        # import subgrid_calculator as sc
        
        # I want to really dig into speed on this code and figure out how to make the
        # the code fast
        
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
            elementDict["DEM%s"%i] = totalEleInfoTable[whichAreInside,0].astype('int')
            # delete those elements fom the total list
            # totalEleInfoTable = np.delete(totalEleInfoTable,whichAreInside,axis=0)
            
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
        
        # now create a surface elevation array
        
        ds = 0.2 # will want to experiment with this and speed
        # surface elevation array for caoluations
        surfaceElevations = np.round(np.arange(-5,5+ds,ds),2).astype('float32') 
        
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
        
        maxElevationEle = np.zeros(numEle).astype(np.float32) # find highest elevation in each element for use in variable phi
        
        maxElevationSubEle = np.zeros((numEle,3)).astype(np.float32) # find the highest elevation in each subElement
        
        # now fill the rows and columns of non subgrid vertices and elements with -99999
        
        wetFraction[np.where(binaryElementList == 0),:,:] = -99999
        area[np.where(binaryElementList == 0)] = -99999
        totWatDepth[np.where(binaryElementList == 0),:,:] = -99999
        wetTotWatDepth[np.where(binaryElementList == 0),:,:] = -99999
        cf[np.where(binaryElementList == 0),:,:] = -99999
        
        # fill max elevation
        maxElevationEle[np.where(binaryElementList == 0)] = -99999
        maxElevationSubEle[np.where(binaryElementList == 0)] = -99999
        
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
                
                # now get quad sub area perimeters
                
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
                    demBathyTopoCut = cp.array(bathyTopo[minCol:maxCol+1,minRow:maxRow+1])
                    manningsnCut = cp.array(manningsn[minCol:maxCol+1,minRow:maxRow+1])
                    
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
                    
                    cupyinsideSubElement = cp.array(insideSubElement)
                    cellsInSubElement = cp.count_nonzero(cupyinsideSubElement)
                    
                    # if there are no cells within the element the DEM is too coarse
                    # you must decrease the DEM resolution in this area
                    if cellsInSubElement == 0:
    
                        sys.exit('DEM {0} resolution too coarse!'.format(i))
                    
                    # get just he bathy topo inside the sub element 
                    
                    bathyTopoInsideSubElement = demBathyTopoCut*cupyinsideSubElement
                    
                    # set 0 values to nan for calculations
                    
                    bathyTopoInsideSubElement[bathyTopoInsideSubElement==0] = cp.nan
                    
                    # get maximum elevation inside the sub element
                    maxElevationSubEle[ele,j] = cp.asnumpy(cp.nanmax(bathyTopoInsideSubElement))
                    
                    # get area of sub element
                    
                    area[ele,j] = tri0Area + tri1Area
                    
                    # print(area[ele,j])
                    
                    # vectorize for surface elevation loop to speed up calcs
                    # remove cells not inside sub element which will flatten the array
                    bathyTopoInsideSubElementNoNaN = bathyTopoInsideSubElement[~cp.isnan(bathyTopoInsideSubElement)]
                    manningsnCutNoNaN = manningsnCut[~cp.isnan(bathyTopoInsideSubElement)]
                    
                    # get the total water depth at each surface elevation
                    
                    temptotWatDepth =  cp.array(surfaceElevations[:,None]) - bathyTopoInsideSubElementNoNaN
                    
                    # count the number of wet cells
                    
                    wetCellsInSubArea = temptotWatDepth > 0.001
                    
                    wetCellsInSubAreaCount = cp.sum(wetCellsInSubArea,axis=1)
                    
                    # now set tot water depth of dry cells to nan
                    
                    temptotWatDepth[temptotWatDepth < 0.001] = cp.nan
                    
                    # add to wet frac array
                    
                    cupywetFraction = wetCellsInSubAreaCount/cellsInSubElement
                                    
                    wetFraction[ele,j,:] = cp.asnumpy(cupywetFraction)
                    
                    # add to total water depth array
                    
                    cupytotWatDepth = cp.nansum(temptotWatDepth,axis=1)/cellsInSubElement
                                    
                    totWatDepth[ele,j,:] = cp.asnumpy(cupytotWatDepth)
                    
                    # get wet total water depth and coefficient of friction
                    
                    # find the mannings for only wet areas then 0 the rest for 
                    # use in calculations 
                    
                    manningsnCutNoNaNWet = manningsnCutNoNaN * wetCellsInSubArea
                    
                    tempcf = (9.81*manningsnCutNoNaNWet**2)/(temptotWatDepth**(1/3))
                    
                    # get wet averaged total water depth
                    
                    cupywetTotWatDepth = cp.nansum(temptotWatDepth,axis=1)/wetCellsInSubAreaCount
                    wetTotWatDepth[ele,j,:] = cp.asnumpy(cupywetTotWatDepth)
                    
                    # get bottom friction
                    cupycf = cp.nansum(tempcf,axis=1)/wetCellsInSubAreaCount # this is correct I need <Cf>W
                    # cupycf = cp.nansum(tempcf,axis=1)/cellsInSubElement # this is incorrect
                    cf[ele,j,:] = cp.asnumpy(cupycf)
                    
                    return cf[ele,j,:],tempcf,manningsnCutNoNaNWet,temptotWatDepth
                    
                    # get rv for advection correction and bottom friction correction
                    cupyrv = cupywetTotWatDepth/(cp.nansum((temptotWatDepth**(3/2))*(tempcf**(-1/2)),axis=1)/wetCellsInSubAreaCount)        
                    rv[ele,j,:] = cp.asnumpy(cupyrv)
                    
                    # get advection correction
                    cupyadv = (1/cupywetTotWatDepth)*(cp.nansum(temptotWatDepth**2/tempcf,axis=1)/wetCellsInSubAreaCount)*cupyrv**2
                    cadv[ele,j,:] = cp.asnumpy(cupyadv)
                    
                    # get corrected bottom friction for level 1 corrections
                    cupycmf = cupywetTotWatDepth*cupyrv**2 # this is corrrect I need <H>W * Rv**2
                    # cupycmf = cupytotWatDepth*cupyrv**2 # this is incorrect
                    cmf[ele,j,:] = cp.asnumpy(cupycmf)
            
                
                # get the maximum elevation inside the element
                maxElevationEle[ele] = np.max(maxElevationSubEle[ele,:])
                
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
        
        # now put sub-elemental quantities on the vertices
        
        wetFractionVertex = np.zeros((numNode,num_SfcElevs))
        wetTotWatDepthVertex = np.zeros((numNode,num_SfcElevs))
        gridTotWatDepthVertex = np.zeros((numNode,num_SfcElevs))
        vertexArea = np.zeros((numNode,num_SfcElevs))
        cfVertex = np.zeros((numNode,num_SfcElevs))
        cmfVertex = np.zeros((numNode,num_SfcElevs))
        # cadvVertex = np.zeros((numNode,num_SfcElevs))
        maxElevationVertex = np.ones(numNode)*-99999
        
        # now fill non subgrid spaces with -99999
        
        wetFractionVertex[np.where(binaryVertexList == 0),:] = -99999
        wetTotWatDepthVertex[np.where(binaryVertexList == 0),:] = -99999
        gridTotWatDepthVertex[np.where(binaryVertexList == 0),:] = -99999
        vertexArea[np.where(binaryVertexList == 0),:] = -99999
        cfVertex[np.where(binaryVertexList == 0),:] = -99999
        maxElevationVertex[np.where(binaryVertexList == 0)] = -99999
        
        # loop through the elements and sum sub-element quantities surrounding each vertex
        
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
                
                # add max elevation
                if(maxElevationVertex[nm0] < maxElevationSubEle[i,0]):
                    
                    maxElevationVertex[nm0] = maxElevationSubEle[i,0]
            
                if(maxElevationVertex[nm1] < maxElevationSubEle[i,1]):
                    
                    maxElevationVertex[nm1] = maxElevationSubEle[i,1]
                    
                if(maxElevationVertex[nm2] < maxElevationSubEle[i,2]):
                    
                    maxElevationVertex[nm2] = maxElevationSubEle[i,2]
                
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
                
                cadv0 = cadv[i,0,:]
                cadv1 = cadv[i,1,:]
                cadv2 = cadv[i,2,:]
            
                cmfVertex[nm0,:] += cmf0 * area[i,0]
                cmfVertex[nm1,:] += cmf1 * area[i,1]
                cmfVertex[nm2,:] += cmf2 * area[i,2]
            
                # cadvVertex[nm0,:] += cadv0 * area[i,0]
                # cadvVertex[nm1,:] += cadv1 * area[i,1]
                # cadvVertex[nm2,:] += cadv2 * area[i,2]
        
        # now average all of these by the vertex areas
        wetFractionVertex[np.where(binaryVertexList == 1)] = wetFractionVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        wetTotWatDepthVertex[np.where(binaryVertexList == 1)] = wetTotWatDepthVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        gridTotWatDepthVertex[np.where(binaryVertexList == 1)] = gridTotWatDepthVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        cfVertex[np.where(binaryVertexList == 1)] = cfVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        cmfVertex[np.where(binaryVertexList == 1)] = cmfVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        # cadvVertex[np.where(binaryVertexList == 1)] = cadvVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        
        # now we need to check to see if there are any vertices in the subgrid
        # but not connected to any elements in the subgrid 
        
        for i in range(numNode):
            
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
                    cmfVertex[i,:] = -99999
                    # cadvVertex[i,:] = -99999
                        
            # otherwise if vertex is outside subgrid area
            # just make everything -99999
            else:
                
                # change the vertex variables to be -99999
                wetFractionVertex[i,:] = -99999
                wetTotWatDepthVertex[i,:] = -99999
                gridTotWatDepthVertex[i,:] = -99999
                cfVertex[i,:] = -99999
                cmfVertex[i,:] = -99999
                # cadvVertex[i,:] = -99999
                
        # now we need to simplify this lookup table
        
        desiredPhiList = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        
        # transpose the dimensions of the elemental subgrid arrays for writing 
        # netCDF
        
        wetFraction = wetFraction.T
        totWatDepth = totWatDepth.T
        cadv = cadv.T
        
        
        # create empty arrays for the reduced tables
        
        depthsEleForLookup = np.zeros((11,3,len(wetFraction[0,0,:])))
        HEleForLookup = np.zeros((11,3,len(wetFraction[0,0,:])))
        cadvForLookup = np.zeros((11,3,len(wetFraction[0,0,:])))
                            
        
        # now loop through the elements 
        
        for i in range(len(wetFraction[0,0,:])):
            
            for j in range(3):
                
                currPhiArray = wetFraction[:,j,i]
                
                # check to see if the element is always fully wet
                
                if(all(currPhiArray == 1.0)):
                    
                    # set all to really low depths so that in the code you know to set at 1
                
                    depthsEleForLookup[:,j,i] = np.min(surfaceElevations)
                    
                    # set the total water depth to the first value so that we have a reference
                    # point to add zeta to
                    
                    HEleForLookup[:,j,i] = totWatDepth[0,j,i]
                    
                    # set cadvection to the last value in the lookup table
                    
                    cadvForLookup[:,j,i] = cadv[-1,j,i]
                    
                # check to see if the element is always fully dry
                    
                elif(all(currPhiArray == 0.0)):
                    
                    # set all to really high dpeths so that in code you know to set to 0.0
                    
                    depthsEleForLookup[:,j,i] = np.max(surfaceElevations) 
                    HEleForLookup[:,j,i] = 0.0
                    
                    # set cadv to the first value of the lookup table
                    
                    cadvForLookup[:,j,i] = cadv[0,j,i]
        
                # if the element is partially wet for all values
                            
                elif(all(currPhiArray > 0.0) and all(currPhiArray < 1.0)):
                    
                    # set to always dry
                    
                    depthsEleForLookup[:,j,i] = np.max(surfaceElevations) 
                    HEleForLookup[:,j,i] = 0.0
                    cadvForLookup[:,j,i] = cadv[0,j,i]
                            
                # if partially wet to fully wet
                
                elif(all(currPhiArray > 0.0) and any(currPhiArray == 1.0)):
                    
                    # set to always wet
                    
                    depthsEleForLookup[:,j,i] = np.min(surfaceElevations)
                    
                    # set the total water depth to the first value so that we have a reference
                    # point to add zeta to
                    
                    HEleForLookup[:,j,i] = totWatDepth[0,j,i]
                    
                    # set cadvection to the last value in the lookup table
                    
                    cadvForLookup[:,j,i] = cadv[-1,j,i]
                            
                # if fully dry to partially wet
                
                elif(all(currPhiArray < 1.0) and any(currPhiArray == 0.0)):
                    
                    # set to always dry
                    
                    depthsEleForLookup[:,j,i] = np.max(surfaceElevations)  
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
                            depthsEleForLookup[k,j,i] = surfaceElevations[equalTo]
                            HEleForLookup[k,j,i] = totWatDepth[equalTo,j,i]
                            cadvForLookup[k,j,i] = cadv[equalTo,j,i]
                            
                        # for phi == 1.0 you want exactly where that is in the currPhiArray
                            
                        elif(desiredPhi == 1.0):
                            
                            equalTo = np.where(currPhiArray == desiredPhi)[0][0]
                            depthsEleForLookup[k,j,i] = surfaceElevations[equalTo]
                            HEleForLookup[k,j,i] = totWatDepth[equalTo,j,i]
                            cadvForLookup[k,j,i] = cadv[equalTo,j,i]
                            
                        # now look for the points in between 0 and 1
                            
                        else:
                            
                            greaterThan = np.where(currPhiArray > desiredPhi)[0][0]
                            lessThan = np.where(currPhiArray < desiredPhi)[0][-1]
                            
                            depthsEleForLookup[k,j,i] = (((desiredPhi - currPhiArray[lessThan])
                                                          /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                          *(surfaceElevations[greaterThan] - surfaceElevations[lessThan])
                                                          + (surfaceElevations[lessThan]))
                            
                            HEleForLookup[k,j,i] = (((desiredPhi - currPhiArray[lessThan])
                                                          /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                          *(totWatDepth[greaterThan,j,i] - totWatDepth[lessThan,j,i])
                                                          + (totWatDepth[lessThan,j,i]))
                            
                            cadvForLookup[k,j,i] = (((desiredPhi - currPhiArray[lessThan])
                                                          /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                          *(cadv[greaterThan,j,i] - cadv[lessThan,j,i])
                                                          + (cadv[lessThan,j,i]))
                            
                            
                                     
            # print('Element {0} completed'.format(i))
            
        # fill the elements that are not contained in the subgrid region with -99999
        depthsEleForLookup[:,:,np.where(binaryElementList==0)[0]] = -99999
        HEleForLookup[:,:,np.where(binaryElementList==0)[0]] = -99999
        cadvForLookup[:,:,np.where(binaryElementList==0)[0]] = -99999
        
        ncFile = nc.Dataset(outputFilename, mode = 'w', format = 'NETCDF4')
        
        # create dimensions
        
        ncFile.createDimension('numEle',numEle) # element dimension
        ncFile.createDimension('numVert',3) # number of vertices per element
        ncFile.createDimension('numPhi',11) # number of possible phi values
        ncFile.createDimension('numSfcElevs',len(surfaceElevations)) # number of surface elevations
        ncFile.createDimension('numNode',numNode) # number of nodes in mesh
        
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
        # vertex coefficient of friction level 1
        cmfVarVertex = ncFile.createVariable('cmfVertex',np.float32,
                                              ('numNode','numSfcElevs'))
        
        # elemental advection correction
        cadvVar = ncFile.createVariable('cadv',np.float32,
                                        ('numPhi','numVert','numEle'))
        
        
        phiSet[:] = desiredPhiList
        wetFractionVarVertex[:,:] = wetFractionVertex
        wetFractionVarDepths[:,:,:] = depthsEleForLookup
        areaVar[:,:] = area
        totWatDepthVar[:,:,:] = HEleForLookup
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
        cmfVarVertex[:,:] = cmfVertex
        cadvVar[:,:,:] = cadvForLookup
        
        ncFile.close()
        
##################### SUBGRID GPU CALULATOR V4 ################################
# Now use the v4 code and add GPU functionality

    def calculateSubgridCorrectionv4GPU(controlFilename):
        
        import sys
        import netCDF4 as nc
        import numpy as np
        import time
        import matplotlib.pyplot as plt
        import cupy as cp
        # use re to read control file in per Shintaro recommendation
        import re
        # import subgrid_calculator as sc
        
        # I want to really dig into speed on this code and figure out how to make the
        # the code fast
        
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
            elementDict["DEM%s"%i] = totalEleInfoTable[whichAreInside,0].astype('int')
            # delete those elements fom the total list
            # totalEleInfoTable = np.delete(totalEleInfoTable,whichAreInside,axis=0)
            
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
        
        # now create a surface elevation array
        
        ds = 0.5 # will want to experiment with this and speed
        # surface elevation array for caoluations
        surfaceElevations = np.round(np.arange(-20,20+ds,ds),2).astype('float32') 
        
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
        
        maxElevationEle = np.zeros(numEle).astype(np.float32) # find highest elevation in each element for use in variable phi
        
        maxElevationSubEle = np.zeros((numEle,3)).astype(np.float32) # find the highest elevation in each subElement
        
        # now fill the rows and columns of non subgrid vertices and elements with -99999
        
        wetFraction[np.where(binaryElementList == 0),:,:] = -99999
        area[np.where(binaryElementList == 0)] = -99999
        totWatDepth[np.where(binaryElementList == 0),:,:] = -99999
        wetTotWatDepth[np.where(binaryElementList == 0),:,:] = -99999
        cf[np.where(binaryElementList == 0),:,:] = -99999
        
        # fill max elevation
        maxElevationEle[np.where(binaryElementList == 0)] = -99999
        maxElevationSubEle[np.where(binaryElementList == 0)] = -99999
        
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
                
                # now get quad sub area perimeters
                
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
                    demBathyTopoCut = cp.array(bathyTopo[minCol:maxCol+1,minRow:maxRow+1])
                    manningsnCut = cp.array(manningsn[minCol:maxCol+1,minRow:maxRow+1])
                    
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
                                            
                    cupyinsideSubElement = cp.array(np.where(insideSubElement.flatten()))
                    cellsInSubElement = cp.count_nonzero(cupyinsideSubElement)
                    
                    if cellsInSubElement == 0:
    
                        sys.exit('DEM {0} resolution too coarse!'.format(i))
                        
                    bathyTopoInsideSubElement = cp.array(demBathyTopoCut.flatten()[cupyinsideSubElement])
                    
                    # get maximum elevation inside the sub element
                    
                    maxElevationSubEle[ele,j] = cp.asnumpy(cp.max(bathyTopoInsideSubElement))
   
                    # get area of sub element
                    
                    area[ele,j] = tri0Area + tri1Area
                    
                    # get the mannings array
                    
                    manningsnCut = cp.array(manningsnCut.flatten()[cupyinsideSubElement])
                    
                    # get total water depth at each elevation
                    
                    temptotWatDepth =  cp.array(surfaceElevations[:,None]) - bathyTopoInsideSubElement
                    
                    # count the number of wet cells
                    
                    wetCellsInSubArea = temptotWatDepth > 0.001
                    
                    wetCellsInSubAreaCount = cp.count_nonzero(wetCellsInSubArea,axis=1)
                    
                    # now set tot water depth of dry cells to nan
                    
                    temptotWatDepth = temptotWatDepth * wetCellsInSubArea
                    
                    # add to wet frac array
                                    
                    cupywetFraction = wetCellsInSubAreaCount/cellsInSubElement
                    
                    wetFraction[ele,j,:] = cp.asnumpy(cupywetFraction)
                    
                    # add to total water depth array
                                    
                    cupytotWatDepth= cp.sum(temptotWatDepth,axis=1)/cellsInSubElement
                    
                    totWatDepth[ele,j,:] = cp.asnumpy(cupytotWatDepth)
                    
                    # get the wet total water depth
                    
                    cupywetTotWatDepth = cp.sum(temptotWatDepth,axis=1)/wetCellsInSubAreaCount
                    wetTotWatDepth[ele,j,:] = cp.asnumpy(cupywetTotWatDepth)
                    
                    # get bottom friction
                    
                    # find the mannings for only wet areas then 0 the rest for 
                    # use in calculations 
                    
                    manningsnCutWet = manningsnCut * wetCellsInSubArea
                    
                    tempcf = (9.81*manningsnCutWet**2)/(temptotWatDepth**(1/3))
                    
                    cupycf = cp.nansum(tempcf,axis=1)/wetCellsInSubAreaCount # this is correct I need <Cf>W
                    # cupycf = cp.nansum(tempcf,axis=1)/cellsInSubElement # this is incorrect
                    cf[ele,j,:] = cp.asnumpy(cupycf)
                    
                    # return cf[ele,j,:],tempcf,manningsnCutWet,temptotWatDepth
                    
                    # get rv for advection correction and bottom friction correction
                    cupyrv = cupywetTotWatDepth/(cp.nansum((temptotWatDepth**(3/2))*(tempcf**(-1/2)),axis=1)/wetCellsInSubAreaCount)        
                    rv[ele,j,:] = cp.asnumpy(cupyrv)
                    
                    # get advection correction
                    cupyadv = (1/cupywetTotWatDepth)*(cp.nansum(temptotWatDepth**2/tempcf,axis=1)/wetCellsInSubAreaCount)*cupyrv**2
                    cadv[ele,j,:] = cp.asnumpy(cupyadv)
                    
                    # get corrected bottom friction for level 1 corrections
                    cupycmf = cupywetTotWatDepth*cupyrv**2 # this is corrrect I need <H>W * Rv**2
                    # cupycmf = cupytotWatDepth*cupyrv**2 # this is incorrect
                    cmf[ele,j,:] = cp.asnumpy(cupycmf)
 
                # get the maximum elevation inside the element
                maxElevationEle[ele] = np.max(maxElevationSubEle[ele,:])
                
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
        
        # now put sub-elemental quantities on the vertices
        
        wetFractionVertex = np.zeros((numNode,num_SfcElevs))
        wetTotWatDepthVertex = np.zeros((numNode,num_SfcElevs))
        gridTotWatDepthVertex = np.zeros((numNode,num_SfcElevs))
        vertexArea = np.zeros((numNode,num_SfcElevs))
        cfVertex = np.zeros((numNode,num_SfcElevs))
        cmfVertex = np.zeros((numNode,num_SfcElevs))
        # cadvVertex = np.zeros((numNode,num_SfcElevs))
        maxElevationVertex = np.ones(numNode)*-99999
        
        # now fill non subgrid spaces with -99999
        
        wetFractionVertex[np.where(binaryVertexList == 0),:] = -99999
        wetTotWatDepthVertex[np.where(binaryVertexList == 0),:] = -99999
        gridTotWatDepthVertex[np.where(binaryVertexList == 0),:] = -99999
        vertexArea[np.where(binaryVertexList == 0),:] = -99999
        cfVertex[np.where(binaryVertexList == 0),:] = -99999
        maxElevationVertex[np.where(binaryVertexList == 0)] = -99999
        
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
                
                # add max elevation
                if(maxElevationVertex[nm0] < maxElevationSubEle[i,0]):
                    
                    maxElevationVertex[nm0] = maxElevationSubEle[i,0]
            
                if(maxElevationVertex[nm1] < maxElevationSubEle[i,1]):
                    
                    maxElevationVertex[nm1] = maxElevationSubEle[i,1]
                    
                if(maxElevationVertex[nm2] < maxElevationSubEle[i,2]):
                    
                    maxElevationVertex[nm2] = maxElevationSubEle[i,2]
                
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
                
                cadv0 = cadv[i,0,:]
                cadv1 = cadv[i,1,:]
                cadv2 = cadv[i,2,:]
            
                cmfVertex[nm0,:] += cmf0 * area[i,0]
                cmfVertex[nm1,:] += cmf1 * area[i,1]
                cmfVertex[nm2,:] += cmf2 * area[i,2]
            
                # cadvVertex[nm0,:] += cadv0 * area[i,0]
                # cadvVertex[nm1,:] += cadv1 * area[i,1]
                # cadvVertex[nm2,:] += cadv2 * area[i,2]
                
        
        # now average all of these by the vertex areas
        wetFractionVertex[np.where(binaryVertexList == 1)] = wetFractionVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        wetTotWatDepthVertex[np.where(binaryVertexList == 1)] = wetTotWatDepthVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        gridTotWatDepthVertex[np.where(binaryVertexList == 1)] = gridTotWatDepthVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        cfVertex[np.where(binaryVertexList == 1)] = cfVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        cmfVertex[np.where(binaryVertexList == 1)] = cmfVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        # cadvVertex[np.where(binaryVertexList == 1)] = cadvVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        
        # now we need to check to see if there are any vertices in the subgrid
        # but not connected to any elements in the subgrid
        end = time.time()
        print('Finished Averaging Elemental Quantities to Vertices in {} s'.format(end-start)) 
        
        start = time.time()
        for i in range(numNode):
            
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
                    cmfVertex[i,:] = -99999
                    # cadvVertex[i,:] = -99999
                        
            # otherwise if vertex is outside subgrid area
            # just make everything -99999
            else:
                
                # change the vertex variables to be -99999
                wetFractionVertex[i,:] = -99999
                wetTotWatDepthVertex[i,:] = -99999
                gridTotWatDepthVertex[i,:] = -99999
                cfVertex[i,:] = -99999
                cmfVertex[i,:] = -99999
                # cadvVertex[i,:] = -99999
        
        end = time.time()
        print('Finished Checking if Vertices were a part of Elements inside Subgrid in {} s'.format(end-start))
        
        # now we need to simplify this lookup table
        
        start = time.time()
        
        desiredPhiList = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        
        # transpose the dimensions of the elemental subgrid arrays for writing 
        # netCDF
        
        wetFraction = wetFraction.T
        totWatDepth = totWatDepth.T
        cadv = cadv.T
        
        
        # create empty arrays for the reduced tables
        
        depthsEleForLookup = np.zeros((11,3,len(wetFraction[0,0,:])))
        HEleForLookup = np.zeros((11,3,len(wetFraction[0,0,:])))
        cadvForLookup = np.zeros((11,3,len(wetFraction[0,0,:])))
                            
        
        # now loop through the elements 
        
        for i in range(len(wetFraction[0,0,:])):
            
            for j in range(3):
                
                currPhiArray = wetFraction[:,j,i]
                
                # check to see if the element is always fully wet
                
                if(all(currPhiArray == 1.0)):
                    
                    # set all to really low depths so that in the code you know to set at 1
                
                    depthsEleForLookup[:,j,i] = np.min(surfaceElevations)
                    
                    # set the total water depth to the first value so that we have a reference
                    # point to add zeta to
                    
                    HEleForLookup[:,j,i] = totWatDepth[0,j,i]
                    
                    # set cadvection to the last value in the lookup table
                    
                    cadvForLookup[:,j,i] = cadv[-1,j,i]
                    
                # check to see if the element is always fully dry
                    
                elif(all(currPhiArray == 0.0)):
                    
                    # set all to really high dpeths so that in code you know to set to 0.0
                    
                    depthsEleForLookup[:,j,i] = np.max(surfaceElevations) 
                    HEleForLookup[:,j,i] = 0.0
                    
                    # set cadv to the first value of the lookup table
                    
                    cadvForLookup[:,j,i] = cadv[0,j,i]
        
                # if the element is partially wet for all values
                            
                elif(all(currPhiArray > 0.0) and all(currPhiArray < 1.0)):
                    
                    # set to always dry
                    
                    depthsEleForLookup[:,j,i] = np.max(surfaceElevations) 
                    HEleForLookup[:,j,i] = 0.0
                    cadvForLookup[:,j,i] = cadv[0,j,i]
                            
                # if partially wet to fully wet
                
                elif(all(currPhiArray > 0.0) and any(currPhiArray == 1.0)):
                    
                    # set to always wet
                    
                    depthsEleForLookup[:,j,i] = np.min(surfaceElevations)
                    
                    # set the total water depth to the first value so that we have a reference
                    # point to add zeta to
                    
                    HEleForLookup[:,j,i] = totWatDepth[0,j,i]
                    
                    # set cadvection to the last value in the lookup table
                    
                    cadvForLookup[:,j,i] = cadv[-1,j,i]
                            
                # if fully dry to partially wet
                
                elif(all(currPhiArray < 1.0) and any(currPhiArray == 0.0)):
                    
                    # set to always dry
                    
                    depthsEleForLookup[:,j,i] = np.max(surfaceElevations)  
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
                            depthsEleForLookup[k,j,i] = surfaceElevations[equalTo]
                            HEleForLookup[k,j,i] = totWatDepth[equalTo,j,i]
                            cadvForLookup[k,j,i] = cadv[equalTo,j,i]
                            
                        # for phi == 1.0 you want exactly where that is in the currPhiArray
                            
                        elif(desiredPhi == 1.0):
                            
                            equalTo = np.where(currPhiArray == desiredPhi)[0][0]
                            depthsEleForLookup[k,j,i] = surfaceElevations[equalTo]
                            HEleForLookup[k,j,i] = totWatDepth[equalTo,j,i]
                            cadvForLookup[k,j,i] = cadv[equalTo,j,i]
                            
                        # now look for the points in between 0 and 1
                            
                        else:
                            
                            greaterThan = np.where(currPhiArray > desiredPhi)[0][0]
                            lessThan = np.where(currPhiArray < desiredPhi)[0][-1]
                            
                            depthsEleForLookup[k,j,i] = (((desiredPhi - currPhiArray[lessThan])
                                                          /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                          *(surfaceElevations[greaterThan] - surfaceElevations[lessThan])
                                                          + (surfaceElevations[lessThan]))
                            
                            HEleForLookup[k,j,i] = (((desiredPhi - currPhiArray[lessThan])
                                                          /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                          *(totWatDepth[greaterThan,j,i] - totWatDepth[lessThan,j,i])
                                                          + (totWatDepth[lessThan,j,i]))
                            
                            cadvForLookup[k,j,i] = (((desiredPhi - currPhiArray[lessThan])
                                                          /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                                          *(cadv[greaterThan,j,i] - cadv[lessThan,j,i])
                                                          + (cadv[lessThan,j,i]))
                            
                            
                                     
            # print('Element {0} completed'.format(i))
            
        # fill the elements that are not contained in the subgrid region with -99999
        depthsEleForLookup[:,:,np.where(binaryElementList==0)[0]] = -99999
        HEleForLookup[:,:,np.where(binaryElementList==0)[0]] = -99999
        cadvForLookup[:,:,np.where(binaryElementList==0)[0]] = -99999
        
        end = time.time()
        print('Finished Reducing Elemental Tables in {} s'.format(end-start))
        
        ncFile = nc.Dataset(outputFilename, mode = 'w', format = 'NETCDF4')
        
        # create dimensions
        
        ncFile.createDimension('numEle',numEle) # element dimension
        ncFile.createDimension('numVert',3) # number of vertices per element
        ncFile.createDimension('numPhi',11) # number of possible phi values
        ncFile.createDimension('numSfcElevs',len(surfaceElevations)) # number of surface elevations
        ncFile.createDimension('numNode',numNode) # number of nodes in mesh
        
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
        # vertex coefficient of friction level 1
        cmfVarVertex = ncFile.createVariable('cmfVertex',np.float32,
                                              ('numNode','numSfcElevs'))
        
        # elemental advection correction
        cadvVar = ncFile.createVariable('cadv',np.float32,
                                        ('numPhi','numVert','numEle'))
        
        
        phiSet[:] = desiredPhiList
        wetFractionVarVertex[:,:] = wetFractionVertex
        wetFractionVarDepths[:,:,:] = depthsEleForLookup
        areaVar[:,:] = area
        totWatDepthVar[:,:,:] = HEleForLookup
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
        cmfVarVertex[:,:] = cmfVertex
        cadvVar[:,:,:] = cadvForLookup
        
        ncFile.close()
        
##################### SUBGRID CPU CALULATOR DEV ################################

##################### SUBGRID GPU CALULATOR DEV ################################
# Now use the v4 code and add GPU functionality

    def calculateSubgridCorrectionDevGPU(controlFilename):
        
        import sys
        import netCDF4 as nc
        import numpy as np
        import time
        import matplotlib.pyplot as plt
        import cupy as cp
        # use re to read control file in per Shintaro recommendation
        import re
        # import subgrid_calculator as sc
        
        # I want to really dig into speed on this code and figure out how to make the
        # the code fast
        
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
            elementDict["DEM%s"%i] = totalEleInfoTable[whichAreInside,0].astype('int')
            # delete those elements fom the total list
            # totalEleInfoTable = np.delete(totalEleInfoTable,whichAreInside,axis=0)
            
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
        
        # now create a surface elevation array
        
        ds = 0.5 # will want to experiment with this and speed
        # surface elevation array for caoluations
        surfaceElevations = np.round(np.arange(-20,20+ds,ds),2).astype('float32') 
        
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
        
        maxElevationEle = np.zeros(numEle).astype(np.float32) # find highest elevation in each element for use in variable phi
        
        maxElevationSubEle = np.zeros((numEle,3)).astype(np.float32) # find the highest elevation in each subElement
        
        # now fill the rows and columns of non subgrid vertices and elements with -99999
        
        wetFraction[np.where(binaryElementList == 0),:,:] = -99999
        area[np.where(binaryElementList == 0)] = -99999
        totWatDepth[np.where(binaryElementList == 0),:,:] = -99999
        wetTotWatDepth[np.where(binaryElementList == 0),:,:] = -99999
        cf[np.where(binaryElementList == 0),:,:] = -99999
        
        # fill max elevation
        maxElevationEle[np.where(binaryElementList == 0)] = -99999
        maxElevationSubEle[np.where(binaryElementList == 0)] = -99999
        
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
                
                # now get quad sub area perimeters
                
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
                    demBathyTopoCut = cp.array(bathyTopo[minCol:maxCol+1,minRow:maxRow+1])
                    manningsnCut = cp.array(manningsn[minCol:maxCol+1,minRow:maxRow+1])
                    
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
                                            
                    cupyinsideSubElement = cp.array(np.where(insideSubElement.flatten()))
                    cellsInSubElement = cp.count_nonzero(cupyinsideSubElement)
                    
                    if cellsInSubElement == 0:
    
                        sys.exit('DEM {0} resolution too coarse!'.format(i))
                        
                    bathyTopoInsideSubElement = cp.array(demBathyTopoCut.flatten()[cupyinsideSubElement])
                    
                    # get maximum elevation inside the sub element
                    
                    maxElevationSubEle[ele,j] = cp.asnumpy(cp.max(bathyTopoInsideSubElement))
   
                    # get area of sub element
                    
                    area[ele,j] = tri0Area + tri1Area
                    
                    # get the mannings array
                    
                    manningsnCut = cp.array(manningsnCut.flatten()[cupyinsideSubElement])
                    
                    # get total water depth at each elevation
                    
                    temptotWatDepth =  cp.array(surfaceElevations[:,None]) - bathyTopoInsideSubElement
                    
                    # count the number of wet cells
                    
                    wetCellsInSubArea = temptotWatDepth > 0.001
                    
                    wetCellsInSubAreaCount = cp.count_nonzero(wetCellsInSubArea,axis=1)
                    
                    # now set tot water depth of dry cells to nan
                    
                    temptotWatDepth = temptotWatDepth * wetCellsInSubArea
                    
                    # add to wet frac array
                                    
                    cupywetFraction = wetCellsInSubAreaCount/cellsInSubElement
                    
                    wetFraction[ele,j,:] = cp.asnumpy(cupywetFraction)
                    
                    # add to total water depth array
                                    
                    cupytotWatDepth= cp.sum(temptotWatDepth,axis=1)/cellsInSubElement
                    
                    totWatDepth[ele,j,:] = cp.asnumpy(cupytotWatDepth)
                    
                    # get the wet total water depth
                    
                    cupywetTotWatDepth = cp.sum(temptotWatDepth,axis=1)/wetCellsInSubAreaCount
                    wetTotWatDepth[ele,j,:] = cp.asnumpy(cupywetTotWatDepth)
                    
                    # get bottom friction
                    
                    # find the mannings for only wet areas then 0 the rest for 
                    # use in calculations 
                    
                    manningsnCutWet = manningsnCut * wetCellsInSubArea
                    
                    tempcf = (9.81*manningsnCutWet**2)/(temptotWatDepth**(1/3))
                    
                    cupycf = cp.nansum(tempcf,axis=1)/wetCellsInSubAreaCount # this is correct I need <Cf>W
                    # cupycf = cp.nansum(tempcf,axis=1)/cellsInSubElement # this is incorrect
                    cf[ele,j,:] = cp.asnumpy(cupycf)
                    
                    # return cf[ele,j,:],tempcf,manningsnCutWet,temptotWatDepth
                    
                    # get rv for advection correction and bottom friction correction
                    cupyrv = cupywetTotWatDepth/(cp.nansum((temptotWatDepth**(3/2))*(tempcf**(-1/2)),axis=1)/wetCellsInSubAreaCount)        
                    rv[ele,j,:] = cp.asnumpy(cupyrv)
                    
                    # get advection correction
                    cupyadv = (1/cupywetTotWatDepth)*(cp.nansum(temptotWatDepth**2/tempcf,axis=1)/wetCellsInSubAreaCount)*cupyrv**2
                    cadv[ele,j,:] = cp.asnumpy(cupyadv)
                    
                    # get corrected bottom friction for level 1 corrections
                    cupycmf = cupywetTotWatDepth*cupyrv**2 # this is corrrect I need <H>W * Rv**2
                    # cupycmf = cupytotWatDepth*cupyrv**2 # this is incorrect
                    cmf[ele,j,:] = cp.asnumpy(cupycmf)
 
                # get the maximum elevation inside the element
                maxElevationEle[ele] = np.max(maxElevationSubEle[ele,:])
                
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
        
        # now put sub-elemental quantities on the vertices
        
        wetFractionVertex = np.zeros((numNode,num_SfcElevs))
        wetTotWatDepthVertex = np.zeros((numNode,num_SfcElevs))
        gridTotWatDepthVertex = np.zeros((numNode,num_SfcElevs))
        vertexArea = np.zeros((numNode,num_SfcElevs))
        cfVertex = np.zeros((numNode,num_SfcElevs))
        cmfVertex = np.zeros((numNode,num_SfcElevs))
        # cadvVertex = np.zeros((numNode,num_SfcElevs))
        maxElevationVertex = np.ones(numNode)*-99999
        
        # now fill non subgrid spaces with -99999
        
        wetFractionVertex[np.where(binaryVertexList == 0),:] = -99999
        wetTotWatDepthVertex[np.where(binaryVertexList == 0),:] = -99999
        gridTotWatDepthVertex[np.where(binaryVertexList == 0),:] = -99999
        vertexArea[np.where(binaryVertexList == 0),:] = -99999
        cfVertex[np.where(binaryVertexList == 0),:] = -99999
        maxElevationVertex[np.where(binaryVertexList == 0)] = -99999
        
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
                
                # add max elevation
                if(maxElevationVertex[nm0] < maxElevationSubEle[i,0]):
                    
                    maxElevationVertex[nm0] = maxElevationSubEle[i,0]
            
                if(maxElevationVertex[nm1] < maxElevationSubEle[i,1]):
                    
                    maxElevationVertex[nm1] = maxElevationSubEle[i,1]
                    
                if(maxElevationVertex[nm2] < maxElevationSubEle[i,2]):
                    
                    maxElevationVertex[nm2] = maxElevationSubEle[i,2]
                
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
                
                # cadv0 = cadv[i,0,:]
                # cadv1 = cadv[i,1,:]
                # cadv2 = cadv[i,2,:]
            
                cmfVertex[nm0,:] += cmf0 * area[i,0]
                cmfVertex[nm1,:] += cmf1 * area[i,1]
                cmfVertex[nm2,:] += cmf2 * area[i,2]
            
                # cadvVertex[nm0,:] += cadv0 * area[i,0]
                # cadvVertex[nm1,:] += cadv1 * area[i,1]
                # cadvVertex[nm2,:] += cadv2 * area[i,2]
        
        # now average all of these by the vertex areas
        wetFractionVertex[np.where(binaryVertexList == 1)] = wetFractionVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        wetTotWatDepthVertex[np.where(binaryVertexList == 1)] = wetTotWatDepthVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        gridTotWatDepthVertex[np.where(binaryVertexList == 1)] = gridTotWatDepthVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        cfVertex[np.where(binaryVertexList == 1)] = cfVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        cmfVertex[np.where(binaryVertexList == 1)] = cmfVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        # cadvVertex[np.where(binaryVertexList == 1)] = cadvVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        end = time.time()
        print('Put elemental quantities to vertices {} s'.format(end-start))
        # now we need to check to see if there are any vertices in the subgrid
        # but not connected to any elements in the subgrid 
        start = time.time()
        
        meshConnectivityInSubgrid = meshConnectivity[containedElementList0Index]
        
        for vertex in containedVertexList0Index:
            
            inSubgrid = np.any(meshConnectivityInSubgrid == vertex)
            
            if(inSubgrid != True):
                # change the vertex to not be in the subgrid
                binaryVertexList[vertex] = 0
                # change the vertex variables to be -99999
                wetFractionVertex[vertex,:] = -99999
                wetTotWatDepthVertex[vertex,:] = -99999
                gridTotWatDepthVertex[vertex,:] = -99999
                cfVertex[vertex,:] = -99999
                cmfVertex[vertex,:] = -99999
                
                
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
        
        
        # now loop through the elements
        
        # first fill the elements that are either all dry or all wet and we will
        # interpolate the ones in between
        # go ahead and set any subelement that does not get fully wet to dry
        checkwhere = np.all(wetFraction<1.0,axis=0).nonzero()
        checkwhereEle = checkwhere[1]
        checkwhereVert  = checkwhere[0]
        depthsEleForLookup[:,checkwhereVert,checkwhereEle] = maxSurfElev
        HEleForLookup[:,checkwhereVert,checkwhereEle] = 0.0
        cadvForLookup[:,checkwhereVert,checkwhereEle] = cadv[0,checkwhereVert,checkwhereEle]
        # go ahead and set any element that gets fully wet to wet and we will interpolate later
        checkwhere1 = np.any(wetFraction==1.0,axis=0).nonzero()
        checkwhereEle1 = checkwhere1[1]
        checkwhereVert1  = checkwhere1[0]
        depthsEleForLookup[:,checkwhereVert1,checkwhereEle1] = minSurfElev
        HEleForLookup[:,checkwhereVert1,checkwhereEle1] = totWatDepth[0,checkwhereVert1,checkwhereEle1]
        cadvForLookup[:,checkwhereVert1,checkwhereEle1] = cadv[-1,checkwhereVert1,checkwhereEle1]
        # now find where there are dry to partially wet to fully wet subelement
        checkwhere1 = np.any(wetFraction== 1.0,axis=0).nonzero()
        chechwhereVert1 = checkwhere1[0]
        checkwhereEle1 = checkwhere1[1]
        checkwhere0 = np.any(wetFraction == 0.0,axis=0).nonzero()
        chechwhereVert0 = checkwhere0[0]
        checkwhereEle0 = checkwhere0[1]
        
        # now find where this overlaps
        
        checkwherebothEle = np.intersect1d(checkwhereEle1,checkwhereEle0)

        intersectNodes = []

        for element in checkwherebothEle:
            
            Nodes0th = chechwhereVert0[np.where(checkwhereEle0 == element)]
            Nodes1st = chechwhereVert1[np.where(checkwhereEle1 == element)]
            intersectNodes.append(np.intersect1d(Nodes0th,Nodes1st))
            
            
        for i in range(len(checkwherebothEle)):
            
            element = checkwherebothEle[i]
            verts = intersectNodes[i]
            
            for vert in verts:

                currPhiArray = wetFraction[:,vert,element]
                
                for k in range(len(desiredPhiList)):
                
                    desiredPhi = desiredPhiList[k]
                    
                    # for phi == 0.0 you want exactly where that is in the currPhiArray
                    
                    if(desiredPhi == 0.0):
                        
                        equalTo = np.where(currPhiArray == desiredPhi)[0][-1]
                        depthsEleForLookup[k,vert,element] = surfaceElevations[equalTo]
                        HEleForLookup[k,vert,element] = totWatDepth[equalTo,vert,element]
                        cadvForLookup[k,vert,element] = cadv[equalTo,vert,element]
                        
                        
                    # for phi == 1.0 you want exactly where that is in the currPhiArray
                        
                    elif(desiredPhi == 1.0):
                        
                        equalTo = np.where(currPhiArray == desiredPhi)[0][0]
                        depthsEleForLookup[k,vert,element] = surfaceElevations[equalTo]
                        HEleForLookup[k,vert,element] = totWatDepth[equalTo,vert,element]
                        cadvForLookup[k,vert,element] = cadv[equalTo,vert,element]
                        
                        
                    # now look for the points in between 0 and 1
                        
                    else:
                        
                        greaterThan = np.where(currPhiArray > desiredPhi)[0][0]
                        lessThan = np.where(currPhiArray < desiredPhi)[0][-1]
                        
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
                        
            # print('Finished Element {}'.format(element))
         
        end = time.time()
        print('Reduction Finished and took {} s'.format(end-start))

        # fill the elements that are not contained in the subgrid region with -99999
        depthsEleForLookup[:,:,np.where(binaryElementList==0)[0]] = -99999
        HEleForLookup[:,:,np.where(binaryElementList==0)[0]] = -99999
        cadvForLookup[:,:,np.where(binaryElementList==0)[0]] = -99999
        
        ncFile = nc.Dataset(outputFilename, mode = 'w', format = 'NETCDF4')
        
        # create dimensions
        
        ncFile.createDimension('numEle',numEle) # element dimension
        ncFile.createDimension('numVert',3) # number of vertices per element
        ncFile.createDimension('numPhi',11) # number of possible phi values
        ncFile.createDimension('numSfcElevs',len(surfaceElevations)) # number of surface elevations
        ncFile.createDimension('numNode',numNode) # number of nodes in mesh
        
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
        # vertex coefficient of friction level 1
        cmfVarVertex = ncFile.createVariable('cmfVertex',np.float32,
                                              ('numNode','numSfcElevs'))
        
        # elemental advection correction
        cadvVar = ncFile.createVariable('cadv',np.float32,
                                        ('numPhi','numVert','numEle'))
        
        
        phiSet[:] = desiredPhiList
        wetFractionVarVertex[:,:] = wetFractionVertex
        wetFractionVarDepths[:,:,:] = depthsEleForLookup
        areaVar[:,:] = area
        totWatDepthVar[:,:,:] = HEleForLookup
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
        cmfVarVertex[:,:] = cmfVertex
        cadvVar[:,:,:] = cadvForLookup
        
        ncFile.close()
        
        
##################### SUBGRID CPU CALULATOR DEV ################################

    def calculateSubgridCorrectionDev(controlFilename):
        
        import sys
        import netCDF4 as nc
        import numpy as np
        import time
        import matplotlib.pyplot as plt
        # use re to read control file in per Shintaro recommendation
        import re
        # import subgrid_calculator as sc
        
        # I want to really dig into speed on this code and figure out how to make the
        # the code fast
        
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
        
        # now create a surface elevation array
        
        ds = 0.5 # will want to experiment with this and speed
        # surface elevation array for caoluations
        surfaceElevations = np.round(np.arange(-20,20+ds,ds),2).astype('float32') 
        
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
        
        maxElevationEle = np.zeros(numEle).astype(np.float32) # find highest elevation in each element for use in variable phi
        
        maxElevationSubEle = np.zeros((numEle,3)).astype(np.float32) # find the highest elevation in each subElement
        
        # now fill the rows and columns of non subgrid vertices and elements with -99999
        
        wetFraction[np.where(binaryElementList == 0),:,:] = -99999
        area[np.where(binaryElementList == 0)] = -99999
        totWatDepth[np.where(binaryElementList == 0),:,:] = -99999
        wetTotWatDepth[np.where(binaryElementList == 0),:,:] = -99999
        cf[np.where(binaryElementList == 0),:,:] = -99999
        
        # fill max elevation
        maxElevationEle[np.where(binaryElementList == 0)] = -99999
        maxElevationSubEle[np.where(binaryElementList == 0)] = -99999
        
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
                    
                    # get maximum elevation inside the sub element
                    
                    maxElevationSubEle[ele,j] = np.nanmax(bathyTopoInsideSubElement)
                    
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
                    
                    wetCellsInSubArea = temptotWatDepth > 0.001
                    
                    wetCellsInSubAreaCount = np.sum(wetCellsInSubArea,axis=1)
                    
                    # now set tot water depth of dry cells to nan
                    
                    temptotWatDepth[temptotWatDepth < 0.001] = np.nan
                    
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
                    
                    cmf[ele,j,:] = wetTotWatDepth[ele,j,:]*rv[ele,j,:]**2 # this is corrrect I need <H>W * Rv**2
                    # cmf[ele,j,:] = totWatDepth[ele,j,:]*rv[ele,j,:]**2 # this is incorrect
                    
                
                # get the maximum elevation inside the element
                maxElevationEle[ele] = np.max(maxElevationSubEle[ele,:])
            
                                        
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
        
        # now put sub-elemental quantities on the vertices
        
        wetFractionVertex = np.zeros((numNode,num_SfcElevs))
        wetTotWatDepthVertex = np.zeros((numNode,num_SfcElevs))
        gridTotWatDepthVertex = np.zeros((numNode,num_SfcElevs))
        vertexArea = np.zeros((numNode,num_SfcElevs))
        cfVertex = np.zeros((numNode,num_SfcElevs))
        cmfVertex = np.zeros((numNode,num_SfcElevs))
        # cadvVertex = np.zeros((numNode,num_SfcElevs))
        maxElevationVertex = np.ones(numNode)*-99999
        
        # now fill non subgrid spaces with -99999
        
        wetFractionVertex[np.where(binaryVertexList == 0),:] = -99999
        wetTotWatDepthVertex[np.where(binaryVertexList == 0),:] = -99999
        gridTotWatDepthVertex[np.where(binaryVertexList == 0),:] = -99999
        vertexArea[np.where(binaryVertexList == 0),:] = -99999
        cfVertex[np.where(binaryVertexList == 0),:] = -99999
        maxElevationVertex[np.where(binaryVertexList == 0)] = -99999
        
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
                
                # add max elevation
                if(maxElevationVertex[nm0] < maxElevationSubEle[i,0]):
                    
                    maxElevationVertex[nm0] = maxElevationSubEle[i,0]
            
                if(maxElevationVertex[nm1] < maxElevationSubEle[i,1]):
                    
                    maxElevationVertex[nm1] = maxElevationSubEle[i,1]
                    
                if(maxElevationVertex[nm2] < maxElevationSubEle[i,2]):
                    
                    maxElevationVertex[nm2] = maxElevationSubEle[i,2]
                
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
                
                # cadv0 = cadv[i,0,:]
                # cadv1 = cadv[i,1,:]
                # cadv2 = cadv[i,2,:]
            
                cmfVertex[nm0,:] += cmf0 * area[i,0]
                cmfVertex[nm1,:] += cmf1 * area[i,1]
                cmfVertex[nm2,:] += cmf2 * area[i,2]
            
                # cadvVertex[nm0,:] += cadv0 * area[i,0]
                # cadvVertex[nm1,:] += cadv1 * area[i,1]
                # cadvVertex[nm2,:] += cadv2 * area[i,2]
        
        # now average all of these by the vertex areas
        wetFractionVertex[np.where(binaryVertexList == 1)] = wetFractionVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        wetTotWatDepthVertex[np.where(binaryVertexList == 1)] = wetTotWatDepthVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        gridTotWatDepthVertex[np.where(binaryVertexList == 1)] = gridTotWatDepthVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        cfVertex[np.where(binaryVertexList == 1)] = cfVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        cmfVertex[np.where(binaryVertexList == 1)] = cmfVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
        # cadvVertex[np.where(binaryVertexList == 1)] = cadvVertex[np.where(binaryVertexList == 1)]/vertexArea[np.where(binaryVertexList == 1)]
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
        
        
        # for vertex in containedVertexList0Index:
            
        #     inSubgrid = np.any(meshConnectivityInSubgrid == vertex)
            
        #     if(inSubgrid != True):
        #         # change the vertex to not be in the subgrid
        #         binaryVertexList[vertex] = 0
        #         # change the vertex variables to be -99999
        #         wetFractionVertex[vertex,:] = -99999
        #         wetTotWatDepthVertex[vertex,:] = -99999
        #         gridTotWatDepthVertex[vertex,:] = -99999
        #         cfVertex[vertex,:] = -99999
        #         cmfVertex[vertex,:] = -99999
                
                
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
        
        
        # now loop through the elements
        
        # first fill the elements that are either all dry or all wet and we will
        # interpolate the ones in between
        # go ahead and set any subelement that does not get fully wet to dry
        checkwhere = np.all(wetFraction<1.0,axis=0).nonzero()
        checkwhereEle = checkwhere[1]
        checkwhereVert  = checkwhere[0]
        depthsEleForLookup[:,checkwhereVert,checkwhereEle] = maxSurfElev
        HEleForLookup[:,checkwhereVert,checkwhereEle] = 0.0
        cadvForLookup[:,checkwhereVert,checkwhereEle] = cadv[0,checkwhereVert,checkwhereEle]
        # go ahead and set any element that gets fully wet to wet and we will interpolate later
        checkwhere1 = np.any(wetFraction==1.0,axis=0).nonzero()
        checkwhereEle1 = checkwhere1[1]
        checkwhereVert1  = checkwhere1[0]
        depthsEleForLookup[:,checkwhereVert1,checkwhereEle1] = minSurfElev
        HEleForLookup[:,checkwhereVert1,checkwhereEle1] = totWatDepth[0,checkwhereVert1,checkwhereEle1]
        cadvForLookup[:,checkwhereVert1,checkwhereEle1] = cadv[-1,checkwhereVert1,checkwhereEle1]
        # now find where there are dry to partially wet to fully wet subelement
        # checkwhere1 = np.any(wetFraction== 1.0,axis=0).nonzero()
        # checkwhereVert1 = checkwhere1[0]
        # checkwhereEle1 = checkwhere1[1]
        checkwhere0 = np.any(wetFraction == 0.0,axis=0).nonzero()
        checkwhereVert0 = checkwhere0[0]
        checkwhereEle0 = checkwhere0[1]
        
        # now find where this overlaps
        # return checkwhere0,checkwhere1
        # checkwherebothEle = np.intersect1d(checkwhereEle1,checkwhereEle0)
        end = time.time()
        print('Finished prepping for Reduction took {} s'.format(end-start))
        start = time.time()
        
        for i in range(len(checkwhereEle0)):
            
            element = checkwhereEle0[i]
            vert = checkwhereVert0[i]
            currPhiArray = wetFraction[:,vert,element]
            
            # make sure that the phi array also gets fully wet and then proceed
            # otherwise just skip
            
            if(1.0 in currPhiArray):
                
                # do when phi == 0
                # for phi == 0.0 you want exactly where that is in the currPhiArray
                equalTo = np.where(currPhiArray == 0.0)[0][-1]
                depthsEleForLookup[0,vert,element] = surfaceElevations[equalTo]
                HEleForLookup[0,vert,element] = totWatDepth[equalTo,vert,element]
                cadvForLookup[0,vert,element] = cadv[equalTo,vert,element]
                
                # do when phi == 1
                # for phi == 1.0 you want exactly where that is in the currPhiArray
                equalTo = np.where(currPhiArray == 1.0)[0][0]
                depthsEleForLookup[-1,vert,element] = surfaceElevations[equalTo]
                HEleForLookup[-1,vert,element] = totWatDepth[equalTo,vert,element]
                cadvForLookup[-1,vert,element] = cadv[equalTo,vert,element]
                
                # only iterate between the second and next to last phi
                for k in range(1,len(desiredPhiList)-1):
                
                    desiredPhi = desiredPhiList[k]
                    greaterThan = np.where(currPhiArray > desiredPhi)[0][0]
                    lessThan = np.where(currPhiArray < desiredPhi)[0][-1]
                    
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
                    
            
            # print('Finished Element {} of {}'.format(i,len(checkwhereEle0)))
                    # # for phi == 0.0 you want exactly where that is in the currPhiArray
                    
                    # if(desiredPhi == 0.0):
                        
                    #     equalTo = np.where(currPhiArray == desiredPhi)[0][-1]
                    #     depthsEleForLookup[k,vert,element] = surfaceElevations[equalTo]
                    #     HEleForLookup[k,vert,element] = totWatDepth[equalTo,vert,element]
                    #     cadvForLookup[k,vert,element] = cadv[equalTo,vert,element]
                        
                        
                    # # for phi == 1.0 you want exactly where that is in the currPhiArray
                        
                    # elif(desiredPhi == 1.0):
                        
                    #     equalTo = np.where(currPhiArray == desiredPhi)[0][0]
                    #     depthsEleForLookup[k,vert,element] = surfaceElevations[equalTo]
                    #     HEleForLookup[k,vert,element] = totWatDepth[equalTo,vert,element]
                    #     cadvForLookup[k,vert,element] = cadv[equalTo,vert,element]
                        
                        
                    # # now look for the points in between 0 and 1
                        
                    # else:
                        
                    #     greaterThan = np.where(currPhiArray > desiredPhi)[0][0]
                    #     lessThan = np.where(currPhiArray < desiredPhi)[0][-1]
                        
                    #     depthsEleForLookup[k,vert,element] = (((desiredPhi - currPhiArray[lessThan])
                    #                                   /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                    #                                   *(surfaceElevations[greaterThan] - surfaceElevations[lessThan])
                    #                                   + (surfaceElevations[lessThan]))
                        
                    #     HEleForLookup[k,vert,element] = (((desiredPhi - currPhiArray[lessThan])
                    #                                   /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                    #                                   *(totWatDepth[greaterThan,vert,element] - totWatDepth[lessThan,vert,element])
                    #                                   + (totWatDepth[lessThan,vert,element]))
                        
                    #     cadvForLookup[k,vert,element] = (((desiredPhi - currPhiArray[lessThan])
                    #                                   /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                    #                                   *(cadv[greaterThan,vert,element] - cadv[lessThan,vert,element])
                    #                                   + (cadv[lessThan,vert,element]))
        # intersectNodes = []

        # for element in checkwherebothEle:
            
        #     Nodes0th = chechwhereVert0[np.where(checkwhereEle0 == element)]
        #     Nodes1st = chechwhereVert1[np.where(checkwhereEle1 == element)]
        #     intersectNodes.append(np.intersect1d(Nodes0th,Nodes1st))
            
        # end = time.time()
        # print('Finished finding which subelements were partially wet took {} s'.format(end-start))
        # start = time.time()
        # for i in range(len(checkwherebothEle)):
            
        #     element = checkwherebothEle[i]
        #     verts = intersectNodes[i]
            
        #     for vert in verts:

        #         currPhiArray = wetFraction[:,vert,element]
                
        #         for k in range(len(desiredPhiList)):
                
        #             desiredPhi = desiredPhiList[k]
                    
        #             # for phi == 0.0 you want exactly where that is in the currPhiArray
                    
        #             if(desiredPhi == 0.0):
                        
        #                 equalTo = np.where(currPhiArray == desiredPhi)[0][-1]
        #                 depthsEleForLookup[k,vert,element] = surfaceElevations[equalTo]
        #                 HEleForLookup[k,vert,element] = totWatDepth[equalTo,vert,element]
        #                 cadvForLookup[k,vert,element] = cadv[equalTo,vert,element]
                        
                        
        #             # for phi == 1.0 you want exactly where that is in the currPhiArray
                        
        #             elif(desiredPhi == 1.0):
                        
        #                 equalTo = np.where(currPhiArray == desiredPhi)[0][0]
        #                 depthsEleForLookup[k,vert,element] = surfaceElevations[equalTo]
        #                 HEleForLookup[k,vert,element] = totWatDepth[equalTo,vert,element]
        #                 cadvForLookup[k,vert,element] = cadv[equalTo,vert,element]
                        
                        
        #             # now look for the points in between 0 and 1
                        
        #             else:
                        
        #                 greaterThan = np.where(currPhiArray > desiredPhi)[0][0]
        #                 lessThan = np.where(currPhiArray < desiredPhi)[0][-1]
                        
        #                 depthsEleForLookup[k,vert,element] = (((desiredPhi - currPhiArray[lessThan])
        #                                               /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
        #                                               *(surfaceElevations[greaterThan] - surfaceElevations[lessThan])
        #                                               + (surfaceElevations[lessThan]))
                        
        #                 HEleForLookup[k,vert,element] = (((desiredPhi - currPhiArray[lessThan])
        #                                               /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
        #                                               *(totWatDepth[greaterThan,vert,element] - totWatDepth[lessThan,vert,element])
        #                                               + (totWatDepth[lessThan,vert,element]))
                        
        #                 cadvForLookup[k,vert,element] = (((desiredPhi - currPhiArray[lessThan])
        #                                               /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
        #                                               *(cadv[greaterThan,vert,element] - cadv[lessThan,vert,element])
        #                                               + (cadv[lessThan,vert,element]))
                        
        #     print('Finished Element {} of {}'.format(i,len(checkwherebothEle)))
         
        end = time.time()
        print('Reduction of partially wet elements finished and took {} s'.format(end-start))

        # fill the elements that are not contained in the subgrid region with -99999
        depthsEleForLookup[:,:,np.where(binaryElementList==0)[0]] = -99999
        HEleForLookup[:,:,np.where(binaryElementList==0)[0]] = -99999
        cadvForLookup[:,:,np.where(binaryElementList==0)[0]] = -99999
        
        # deallocate arrays
        totWatDepth = None
        cadv = None
        wetFraction = None
        
        ncFile = nc.Dataset(outputFilename, mode = 'w', format = 'NETCDF4')
        
        # create dimensions
        
        ncFile.createDimension('numEle',numEle) # element dimension
        ncFile.createDimension('numVert',3) # number of vertices per element
        ncFile.createDimension('numPhi',11) # number of possible phi values
        ncFile.createDimension('numSfcElevs',len(surfaceElevations)) # number of surface elevations
        ncFile.createDimension('numNode',numNode) # number of nodes in mesh
        
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
        # vertex coefficient of friction level 1
        cmfVarVertex = ncFile.createVariable('cmfVertex',np.float32,
                                              ('numNode','numSfcElevs'))
        
        # elemental advection correction
        cadvVar = ncFile.createVariable('cadv',np.float32,
                                        ('numPhi','numVert','numEle'))
        
        
        phiSet[:] = desiredPhiList
        wetFractionVarVertex[:,:] = wetFractionVertex
        wetFractionVarDepths[:,:,:] = depthsEleForLookup
        areaVar[:,:] = area
        totWatDepthVar[:,:,:] = HEleForLookup
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
        cmfVarVertex[:,:] = cmfVertex
        cadvVar[:,:,:] = cadvForLookup
        
        ncFile.close()
        
        #             # # now combine the two triangles and find the points inside
                    
                    # insideSubElement = np.logical_or(insideTri0,insideTri1)
                    
                    # # count the number of subgrid cells within the subelement
                        
                    # insideSubElement = np.where(insideSubElement.flatten())
                    # cellsInSubElement = np.count_nonzero(insideSubElement)
                    
                    # if cellsInSubElement == 0:
                    
                    #     sys.exit('DEM {0} resolution too coarse!'.format(i))
                        
                    # bathyTopoInsideSubElement = demBathyTopoCut.flatten()[insideSubElement]
                    
                    # # get maximum elevation inside the sub element
                    
                    # maxElevationSubEle[ele,j] = np.max(bathyTopoInsideSubElement)
                    
                    # # get area of sub element
                    
                    # area[ele,j] = tri0Area + tri1Area
                    
                    # # get the mannings array
                    
                    # manningsnCut = manningsnCut.flatten()[insideSubElement]
                    
                    # # get total water depth at each elevation
                    
                    # temptotWatDepth =  surfaceElevations[:,None] - bathyTopoInsideSubElement
                    
                    # # count the number of wet cells
                    
                    # wetCellsInSubArea = temptotWatDepth > 0.001
                    
                    # wetCellsInSubAreaCount = np.count_nonzero(wetCellsInSubArea,axis=1)
                    
                    # # now set tot water depth of dry cells to nan
                    
                    # temptotWatDepth = temptotWatDepth * wetCellsInSubArea
                    
                    # # add to wet frac array
                                    
                    # wetFraction[ele,j,:] = wetCellsInSubAreaCount/cellsInSubElement
                    
                    # # add to total water depth array
                                    
                    # totWatDepth[ele,j,:] = np.sum(temptotWatDepth,axis=1)/cellsInSubElement
                    
                    # # get the wet total water depth
                    # # tempwetTotWatDepth = np.sum(temptotWatDepth,axis=1)/wetCellsInSubAreaCount
                    # wetTotWatDepth[ele,j,:] = np.sum(temptotWatDepth,axis=1)/wetCellsInSubAreaCount
                    
                    # # get bottom friction
                    
                    # # find the mannings for only wet areas then 0 the rest for 
                    # # use in calculations 
                    
                    # manningsnCutWet = manningsnCut * wetCellsInSubArea
                    
                    # tempcf = (9.81*manningsnCutWet**2)/(temptotWatDepth**(1/3))
                    
                    # cf[ele,j,:] = np.nansum(tempcf,axis=1)/wetCellsInSubAreaCount # this is correct I need <Cf>W
                    # # cupycf = cp.nansum(tempcf,axis=1)/cellsInSubElement # this is incorrect
                    
                    # # return cf[ele,j,:],tempcf,manningsnCutWet,temptotWatDepth
                    
                    # # get rv for advection correction and bottom friction correction
                    # # temprv = tempwetTotWatDepth/(np.nansum((temptotWatDepth**(3/2))*(tempcf**(-1/2)),axis=1)/wetCellsInSubAreaCount)
                    # rv[ele,j,:] = wetTotWatDepth[ele,j,:]/(np.nansum((temptotWatDepth**(3/2))*(tempcf**(-1/2)),axis=1)/wetCellsInSubAreaCount)   
                    
                    # # get advection correction
                    # cadv[ele,j,:] = (1/wetTotWatDepth[ele,j,:])*(np.nansum(temptotWatDepth**2/tempcf,axis=1)/wetCellsInSubAreaCount)*rv[ele,j,:]**2
                    
                    # # get corrected bottom friction for level 1 corrections
                    # cmf[ele,j,:] = wetTotWatDepth[ele,j,:]*rv[ele,j,:]**2 # this is corrrect I need <H>W * Rv**2
                    # # cupycmf = cupytotWatDepth*cupyrv**2 # this is incorrect
                    
        
