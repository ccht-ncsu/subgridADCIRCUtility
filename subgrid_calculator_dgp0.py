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
import geopandas as gpd
from shapely.geometry import Point, Polygon

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

    def __init__(self, subgridControlFilename = None, level0andLevel1 = True, GPU = False):
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
        self.surfElevIncrement = 0.1
        self.maxWatDepthAboveHighestGround = 3.0
        self.watDepthAboveHighestGroundIncrement = 1.0

        self.default_manningsn = 0.03
        self.cf_lower_lim = 0.0025

        self.need_projection = True

        if subgridControlFilename != None:
            self.readSubgridControlFile(subgridControlFilename)

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

    ######################## Function to set control parameters ###################
    
    def setControlParameters(self, outputFilename, meshFilename, numDEMs=0, demFilenameList=[], numLCs=0, landcoverFilenameList=[]):
        self.control.outputFilename = outputFilename
        self.control.meshFilename = meshFilename
        self.control.numDEMs = numDEMs
        self.control.demFilenameList = demFilenameList
        self.control.numLCs = numLCs
        self.control.landcoverFilenameList = landcoverFilenameList
        self.control.loaded = True

    ##################### READ IN A MESH IN FORT.14 FORMAT #######################

    def readMeshFile(self, meshFilename=None):
        
        import pandas as pd
        import matplotlib.tri as mtri
        import numpy as np
        
        if meshFilename == None:
            meshFilename = self.control.meshFilename
        else:
            self.control.meshFilename = meshFilename

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

        self.createEdgeList()

        return

    ############ CALCULATE AREA OF A TRIANGLE #######################################
    def build_gpd_elem_triangles(self):
        polys = [Polygon(zip(self.mesh.coord.Longitude[tri],self.mesh.coord.Latitude[tri])) for tri in self.mesh.tri]
        gdf = gpd.GeoDataFrame(index=list(range(len(polys))), crs='epsg:4326', geometry=polys)
        self.mesh.gpd_elem_triangles = gdf

        return

    ############ CALCULATE AREA OF A TRIANGLE #######################################
    
    def triarea(x1,y1,x2,y2,x3,y3):
        
        area = 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    
        return area
    
    ############ CALCULATE LENGTH OF A LINE SEGMENT #######################################
    
    def seglength(x1,y1,x2,y2):
        import math

        length = math.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))
    
        return length
    
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
 
    ############## CHECK IF A DEM CELL IS ALONG A LINE SEGMENT ##################
    
    def isAlong(x1,y1,x2,y2,x,y,xDEMRes,yDEMRes):
        import math
        import sys
        from subgrid_calculator_dgp0 import SubgridCalculatorDGP0 as scm

        dbuf = 1e-8

        dl = math.sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))
        tx = (x2 - x1)/dl
        ty = (y2 - y1)/dl
        if tx == 0.0:
            jd1 = -1
            jd2 = jd1
        elif tx > 0.0:
            jd1 = -1
            jd2 = 0
        else:
            jd1 = 0
            jd2 = -1
        if ty == 0.0:
            id1 = 0
            id2 = id1
        elif ty > 0.0:
            id1 = 0
            id2 = -1
        else:
            id1 = -1
            id2 = 0

        i1 = np.where(np.logical_and(y1 >= y - 0.5*yDEMRes,y1 <= y + 0.5*yDEMRes))[0][id1]
        j1 = np.where(np.logical_and(x1 >= x - 0.5*xDEMRes,x1 <= x + 0.5*xDEMRes))[1][jd1]
        i2 = np.where(np.logical_and(y2 >= y - 0.5*yDEMRes,y2 <= y + 0.5*yDEMRes))[0][id2]
        j2 = np.where(np.logical_and(x2 >= x - 0.5*xDEMRes,x2 <= x + 0.5*xDEMRes))[1][jd2]

        mask = np.zeros_like(x,bool)
        mask[i1,j1] = True

        ni = x.shape[0]
        nj = x.shape[1]

        ii1 = i1
        jj1 = j1

        p0x = x1
        p0y = y1

        while True:
            if tx >= 0.0:
                jj2 = min(jj1 + 1,nj-1)
            else:
                jj2 = max(jj1 - 1,0)
            if ty >= 0.0:
                ii2 = max(ii1 - 1,0)
            else:
                ii2 = min(ii1 + 1,ni-1)

            if ii1 < 0 or ii1 >= ni:
                raise Exception('Invalid i1 value')
            if jj1 < 0 or jj1 >= nj:
                raise Exception('Invalid j1 value')
                   
            while True:
                p1x = x[ii1,jj1]
                p1y = y[ii1,jj1]            

                # Test the upper or lower neighbor cell
                if ii1 != ii2 and ty != 0.0:
                    p2x = x[ii2,jj1]
                    p2y = y[ii2,jj1]
                    p3x = 0.5*(p1x+p2x)
                    p3y = 0.5*(p1y+p2y)                        

                    d0 = abs(p3y - p0y)
                    alpha = d0/abs(ty)

                    p4x = p0x + alpha*tx
                    p4y = p0y + alpha*ty

                    if p4x >= p3x - 0.5*xDEMRes - dbuf and p4x <= p3x + 0.5*xDEMRes +dbuf:
                        mask[ii2,jj1] = True
                        p0x = p4x
                        p0y = p4y
                        ii1 = ii2
                        break

                # Test the left or right neighbor cell
                if jj1 != jj2 and tx != 0.0:
                    p2x = x[ii1,jj2]
                    p2y = y[ii1,jj2]
                    p3x = 0.5*(p1x+p2x)
                    p3y = 0.5*(p1y+p2y)                        

                    d0 = abs(p3x - p0x)
                    alpha = d0/abs(tx)

                    p4x = p0x + alpha*tx
                    p4y = p0y + alpha*ty

                    if p4y >= p3y - 0.5*yDEMRes - dbuf and p4y <= p3y + 0.5*yDEMRes + dbuf:
                        mask[ii1,jj2] = True
                        p0x = p4x
                        p0y = p4y
                        jj1 = jj2
                        break

                raise Exception('Next cell not found')
            
            if ii1 == i2 and jj1 == j2:
                break

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
                    
    ########## CREATE EDGE LIST ##########

    def createEdgeList(self):
        import numpy as np

        numNod = self.mesh.numNod
        numEle = self.mesh.numEle

        # Create the neighbor table
        nndel = np.zeros(numNod,dtype=int)
        for iel in range(numEle):
            nndel[self.mesh.tri[iel,0]] += 1
            nndel[self.mesh.tri[iel,1]] += 1
            nndel[self.mesh.tri[iel,2]] += 1

        mnndel = 0
        for i in range(numNod):
            if(mnndel < nndel[i]):
                mnndel = nndel[i]

        nndel = np.zeros(numNod,dtype=int)
        ndel = np.zeros((numNod,mnndel),dtype=int)
        for iel in range(numEle):
            n1 = self.mesh.tri[iel,0]
            n2 = self.mesh.tri[iel,1]
            n3 = self.mesh.tri[iel,2]

            ndel[n1,nndel[n1]] = iel
            nndel[n1] += 1
            ndel[n2,nndel[n2]] = iel
            nndel[n2] += 1
            ndel[n3,nndel[n3]] = iel
            nndel[n3] += 1

        # Create the edge list
        edflg = np.zeros((numEle,3),dtype=int)
        neled = np.zeros((numEle,3),dtype=int)
        led = np.zeros((3,2),dtype=int)
        nedges = 0

        # Counting the number of edges
        for iel in range(numEle):
            n1 = self.mesh.tri[iel,0]
            n2 = self.mesh.tri[iel,1]
            n3 = self.mesh.tri[iel,2]
            led[0,0] = n2
            led[0,1] = n3
            led[1,0] = n3
            led[1,1] = n1
            led[2,0] = n1
            led[2,1] = n2

            for ied in range(3):
                if(edflg[iel,ied] == 1):
                    continue

                i1 = led[ied,0]
                i2 = led[ied,1]

                nedges = nedges + 1

                for jjel in range(nndel[i1]):
                    jel = ndel[i1,jjel]
                    if jel == iel:
                        continue
                    for jed in range(3):
                        j1 = self.mesh.tri[jel,(jed+1)%3]
                        j2 = self.mesh.tri[jel,(jed+2)%3]
                        if (j1 == i1 and j2 == i2) or (j1 == i2 and j2 == i1):
                            edflg[jel,jed] = 1

        edflg = np.zeros((numEle,3),dtype=int)
        neled = np.zeros((numEle,3),dtype=int)
        nedno = np.zeros((nedges,2),dtype=int)
        nedel = np.full((nedges,2),-1,dtype=int)
        nedges = 0
        for iel in range(numEle):
            n1 = self.mesh.tri[iel,0]
            n2 = self.mesh.tri[iel,1]
            n3 = self.mesh.tri[iel,2]
            led[0,0] = n2
            led[0,1] = n3
            led[1,0] = n3
            led[1,1] = n1
            led[2,0] = n1
            led[2,1] = n2

            for ied in range(3):
                if(edflg[iel,ied] == 1):
                    continue

                i1 = led[ied,0]
                i2 = led[ied,1]

                ed_id = nedges
                nedges = nedges + 1

                neled[iel,ied] = ed_id
                nedno[ed_id,:] = [i1,i2]
                nedel[ed_id,:] = [iel,-1]
                edflg[iel,ied] = 1

                for jjel in range(nndel[i1]):
                    jel = ndel[i1,jjel]
                    if jel == iel:
                        continue
                    for jed in range(3):
                        j1 = self.mesh.tri[jel,(jed+1)%3]
                        j2 = self.mesh.tri[jel,(jed+2)%3]
                        if (j1 == i1 and j2 == i2) or (j1 == i2 and j2 == i1):
                            neled[jel,jed] = ed_id
                            nedel[ed_id,1] = jel
                            edflg[jel,jed] = 1
        
        self.mesh.nnd2el = nndel
        self.mesh.nd2el = ndel
        self.mesh.ed2nd = nedno
        self.mesh.el2ed = neled
        self.mesh.ed2el = nedel
        self.mesh.numEdg = nedges

    ########## CALCULATE ELEMENT SUBGRID CORRECTION FOR VECTORIZED STORAGE ##########

    def calculateElementLookupTableForVectorizedStorage(self):
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

        if self.control.demFilenameList:
            nDEMFilenameList = len(self.control.demFilenameList)
        else:
            nDEMFilenameList = 1

        if self.control.demFilenameList:
            for i in range(nDEMFilenameList):
                
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
        else:
            try:
                # set elevationData from user defined method
                elevationData = self.importDEM_external()
            except:
                raise Exception('Class method importDEM_external not defined') 

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
            # elevationData = None
            
            # get dem bounds to determine what dem to use when looping
            
            elevationDict["bounds%s"%0] = [np.min(xDEMCoordsTemp),
                                            np.max(xDEMCoordsTemp),
                                            np.min(yDEMCoordsTemp),
                                            np.max(yDEMCoordsTemp)]                
             
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
        
        # find if each element is within any of the given polygons
        if hasattr(self.control, 'sgs_region_mask') and self.control.sgs_region_mask:
            polys = gpd.read_file(self.control.sgs_region_mask).unary_union
            elemInsidePolygon = \
            [all(gpd.GeoSeries([
                Point(xS[j,0], yS[j,0]),
                Point(xS[j,1], yS[j,1]),
                Point(xS[j,2], yS[j,2])]
                ).set_crs(epsg="4326").within(polys)) for j in range(numEle)]
        else:
            elemInsidePolygon = np.ones(numEle).astype(bool)

        # loop through DEMs and determine which elements are in them
        # make sure to have DEMs in priority order meaning if you want to use fine
        # resolution in some areas have those dems listed first
        
        containedElementList0Index = []
        noElementDEMs = []
        
        for i in range(nDEMFilenameList):

            # find which elements are in the DEM
            totalEleInfoTable[:,5] = ((totalEleInfoTable[:,1]>=elevationDict["bounds%s"%i][0])
                                      & ((totalEleInfoTable[:,2])<=elevationDict["bounds%s"%i][1])
                                      & ((totalEleInfoTable[:,3])>=elevationDict["bounds%s"%i][2])
                                      & ((totalEleInfoTable[:,4])<=elevationDict["bounds%s"%i][3]))
        
            whichAreInside = list(np.where([totalEleInfoTable[j,5] == 1 and elemInsidePolygon[j] for j in range(numEle)])[0])
            elementDict["DEM%s"%i] = totalEleInfoTable[whichAreInside,0].astype(int)        # store element numbers of the elements inside the DEM bound
            
            whichAreInsideActualEleNumber = totalEleInfoTable[whichAreInside,0].astype(int) # get the actual element numbers 
            totalEleInfoTable = np.delete(totalEleInfoTable,whichAreInside,axis=0)          # delete elements so we can skip those in the next dem

            # keep track if a dem does not have any elements inside to throw and exception later
            if len(whichAreInside) == 0 and self.control.demFilenameList:
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

        # fill area array
        area[np.where(binaryElementList == 0)] = -99999
        
        # fill min/max elevation
        minElevationEle[np.where(binaryElementList == 0)] = 99999
        maxElevationEle[np.where(binaryElementList == 0)] = -99999

        # allocate the index that maps element to location in vector
        elemIndex = np.zeros((numEle)).astype(int)
        elemNumLevel = np.zeros((numEle)).astype(int)
        elemIndexCnt = 0

        # these variables are used if you want level 1 corrections
        if self.level0andLevel1:
            rv = np.array([],dtype=np.float32)
            cmf = np.array([],dtype=np.float32)

        # specify buffer to use in area calculator
        areaDif = 0.00001 
        # areaDif = 0.00000001 
        
        # create variable to keep track of what DEM you have read in
        
        # dictionary to translate between C-CAP and Manning's values
        # landCoverToManning = {0:0.02,2:0.12,3:0.12,4:0.12,5:0.035,6:0.1,7:0.05,8:0.035,
        #                   9:0.16,10:0.18,11:0.17,12:0.08,13:0.15,14:0.075,
        #                   15:0.06,16:0.15,17:0.07,18:0.05,19:0.03,20:0.03,
        #                   21:0.025,22:0.035,23:0.03,25:0.012}
        # change mannings conversion to match OM2D
        # landCoverToManning = {0:0.02, 2:0.15, 3:0.10, 4:0.05, 5:0.02,
        #                         6:0.037, 7:0.033, 8:0.034, 9:0.1, 10:0.11,
        #                         11:0.1, 12:0.05, 13:0.1, 14:0.048, 15:0.045,
        #                         16:0.1, 17:0.048, 18:0.045, 19:0.04,
        #                         20:0.09, 21:0.02, 22:0.015, 23:0.015, 
        #                         24:0.09, 25:0.01}
        # updated to SACS C-CAP mannings table
        # landCoverToManning = {0:0.025, 2:0.12, 3:0.10, 4:0.07, 5:0.035,
        #                       6:0.01, 7:0.055, 8:0.035, 9:0.16, 10:0.18,
        #                       11:0.17, 12:0.08, 13:0.15, 14:0.075, 15:0.07,
        #                       16:0.15, 17:0.07, 18:0.05, 19:0.03,
        #                       20:0.03, 21:0.025, 22:0.035, 23:0.03, 
        #                       24:0.09, 25:0.01}
        
        # landCoverValues = [0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25]
        # add landcover 24
        # landCoverValues = [0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25]
        
        # dictionary to translate between NLCD and Manning's values
        landCoverToManning = {0:0.02, 11:0.02, 12:0.01, 21:0.02, 22:0.05, 23:0.1,
                              24:0.15, 31:0.09, 32:0.04, 41:0.1, 42:0.11, 43:0.1,
                              51:0.04, 52:0.05, 71:0.034, 72:0.03, 73:0.027, 74:0.025,
                              81:0.033, 82:0.037, 90:0.1, 91:0.1, 92:0.048, 93:0.1,
                              94:0.048, 95:0.045, 96:0.045, 97:0.045, 98:0.015,
                              99:0.015, 127:0.02}
        landCoverValues = landCoverToManning.keys()

        # first create a loop for DEMs
        for i in range(nDEMFilenameList):
            
            # reading in DEM again
            # all variables the same as before
            if self.control.demFilenameList:
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

                # convert values
                for value in landCoverValues:
                    nArray[nArray == value] = landCoverToManning[value]
                
                
                # nArray[np.isnan(nArray)] = 0.012  # set nan values to 0.012 for ccap
                
                nArray[np.isnan(nArray)] = 0.02     # set nan values to 0.02
            else:
                if not 'elevationData' in locals():
                    try:
                        # set elevationData from user defined method
                        elevationData = self.importDEM_external()
                    except:
                        raise Exception('Class method importDEM_external not defined') 
                
                xDEM = elevationData[0]
                yDEM = elevationData[1]
                zDEM = elevationData[2]
                xDEMRes = elevationData[3]
                yDEMRes = -1*elevationData[4]
                xDEMCoords = elevationData[5]
                yDEMCoords = elevationData[6]
                # elevationData = None # deallocate 
                
                print("Default Manning's N value = {} is used.".format(self.default_manningsn))
                nArray = np.full(zDEM.shape,self.default_manningsn,dtype=np.double) # array of mannings n values
                            
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
                if self.need_projection:
                    xyCurrElementMeters = scm.projectMeshToMercator(xyArray[:,1],
                                                                    xyArray[:,0])
                else:
                    xyCurrElementMeters = [xyArray[:,0],xyArray[:,1]]
                
                # calculate the area of the element
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
                # print(np.max(currYPerimeterPoints),np.min(currYPerimeterPoints),np.min(currXPerimeterPoints),np.max(currXPerimeterPoints))
                # print(maxY,minY,minX,maxX)
                # print(minRow,maxRow,minCol,maxCol)
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
                mask = np.logical_and(mask, np.isfinite(zCutGeoTiffMatrix2))
                # print(xCutGeoTiffMatrix2)
                # print(yCutGeoTiffMatrix2)
                
                # convert mask to cupy array
                mask = cp.asarray(mask)
                zCutGeoTiffMatrix2 = cp.asarray(zCutGeoTiffMatrix2)
                nCutGeoTiffMatrix2 = cp.asarray(nCutGeoTiffMatrix2)
                
                zCutGeoTiffMatrix2masked = zCutGeoTiffMatrix2[mask]
                nCutGeoTiffMatrix2masked = nCutGeoTiffMatrix2[mask]
                
                # get the min/max elevation of the element
                # print(zCutGeoTiffMatrix2masked.shape,nCutGeoTiffMatrix2masked.shape)
                minElev = cp.nanmin(zCutGeoTiffMatrix2masked)
                maxElev = cp.nanmax(zCutGeoTiffMatrix2masked)

                # count how many cells are within element
                countIn = cp.count_nonzero(mask)
                
                # if there are no cells within the element the DEM is too coarse
                # you must decrease the DEM resolution in this area
                if countIn == 0:
                    sys.exit('DEM {0} resolution too coarse!'.format(i))
            
                # keep track of this for use later
                countInElement += countIn
            
                # create array of surface elevations
                # this is used to calculate the subgrid variables for varying water depths
                surfElevIncrement = self.surfElevIncrement
                minElevFloor = math.floor(minElev/surfElevIncrement)*surfElevIncrement
                surfaceElevations = np.arange(minElevFloor,maxElev+surfElevIncrement,surfElevIncrement)
                if len(surfaceElevations) < 2:
                    surfaceElevations = np.array([minElev,maxElev]).astype(np.float32)
                else:
                    surfaceElevations[0] = minElev
                    surfaceElevations[-1] = maxElev
                r1 = maxElev + self.watDepthAboveHighestGroundIncrement
                r2 = r1 + self.maxWatDepthAboveHighestGround
                surfaceElevations = np.append(surfaceElevations, np.arange(r1,r2))
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
                tempcf[cp.where(tempcf==0.0)] = np.nan
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
                if self.level0andLevel1:
                    rvBottomTermList[cp.where(rvBottomTermList==0.0)] = np.nan
                    cmfTemp = (wetAvgTotWatDepth)*(wetAvgTotWatDepth/(rvBottomTermList/wetDryElementList))**2
                    cmfTemp[cp.where(cmfTemp==np.nan)] = 0.0

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
        cf[cf<self.cf_lower_lim] = self.cf_lower_lim
        cf[np.isnan(cf)] = self.cf_lower_lim
        
        if self.level0andLevel1:
            cmf[cmf<self.cf_lower_lim] = self.cf_lower_lim
            cmf[np.isnan(cmf)] = self.cf_lower_lim
            
        wetFraction[np.isnan(wetFraction)] = 0.0
        wetTotWatDepth[np.isnan(wetTotWatDepth)] = 0.0
        totWatDepth[np.isnan(totWatDepth)] = 0.0

        # sort the vectors in the ascending order of the elemental no.
        print('Sorting...',end='')
        def sortarray(vals):
            sorted = np.zeros(len(vals)).astype(np.float32)
            sortedIndexStart = np.zeros(numEle).astype(int)
            sortedIndexEnd = np.zeros(numEle).astype(int)
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
        self.subgridvectorized.cf = cf
        if self.level0andLevel1:
            self.subgridvectorized.cmf = cmf
        self.subgridvectorized.loaded = True
        
    ########## CALCULATE EDGE SUBGRID CORRECTION FOR VECTORIZED STORAGE ##########

    def calculateEdgeLookupTableForVectorizedStorage(self):
        import sys
        import math
        import numpy as np
        import numpy as cp
        import time
        import matplotlib.pyplot as plt
        from subgrid_calculator_dgp0 import SubgridCalculatorDGP0 as scm

        startTot = time.perf_counter()
        
        # ignore any true divide errors
        np.seterr(invalid='ignore')
        
        # read in mesh
        numNod = self.mesh.numNod
        numEle = self.mesh.numEle
        numEdg = self.mesh.numEdg

        # put edge node coordinates in arrays
        xS = np.vstack((np.asarray(self.mesh.coord['Longitude'])[self.mesh.ed2nd[:,0]],
                        np.asarray(self.mesh.coord['Longitude'])[self.mesh.ed2nd[:,1]])).T
        yS = np.vstack((np.asarray(self.mesh.coord['Latitude'])[self.mesh.ed2nd[:,0]],
                        np.asarray(self.mesh.coord['Latitude'])[self.mesh.ed2nd[:,1]])).T
       
        # create a dictionary to hold the dem data
        elevationDict = {}
        
        if self.control.demFilenameList:
            nDEMFilenameList = len(self.control.demFilenameList)
        else:
            nDEMFilenameList = 1

        if self.control.demFilenameList:
            for i in range(len(self.control.demFilenameList)):            
            
                elevationData = scm.importDEM(self.control.demFilenameList[i])  # read in DEM        
                xDEMTemp = elevationData[0]        # x coordinates of DEM
                yDEMTemp = elevationData[1]        # y coordinates of DEM
                zDEMTemp = elevationData[2]        # elevations of DEM
                xDEMResTemp = elevationData[3]     # resolution in the x direction
                yDEMResTemp = -1*elevationData[4]  # resolution in the y direction
                xDEMCoordsTemp = elevationData[5]  # x coordinates of DEM in 1D array
                yDEMCoordsTemp = elevationData[6]  # y coordinates of DEM in 1D array
                elevationData = None               # deallocate DEM data
                
                # get dem bounds to determine what dem to use when looping            
                elevationDict["bounds%s"%i] = [np.min(xDEMCoordsTemp),
                                            np.max(xDEMCoordsTemp),
                                            np.min(yDEMCoordsTemp),
                                            np.max(yDEMCoordsTemp)]
                
                print('Finished reading DEM {0}.'.format(i))
        else:
            try:
                # set elevationData from user defined method
                elevationData = self.importDEM_external()
            except:
                raise Exception('Class method importDEM_external not defined') 
            
            xDEMTemp = elevationData[0]        # x coordinates of DEM
            yDEMTemp = elevationData[1]        # y coordinates of DEM
            zDEMTemp = elevationData[2]        # elevations of DEM
            xDEMResTemp = elevationData[3]     # resolution in the x direction
            yDEMResTemp = -1*elevationData[4]  # resolution in the y direction
            xDEMCoordsTemp = elevationData[5]  # x coordinates of DEM in 1D array
            yDEMCoordsTemp = elevationData[6]  # y coordinates of DEM in 1D array
            # elevationData = None               # deallocate DEM data
            
            # get dem bounds to determine what dem to use when looping            
            elevationDict["bounds%s"%0] = [np.min(xDEMCoordsTemp),
                                        np.max(xDEMCoordsTemp),
                                        np.min(yDEMCoordsTemp),
                                        np.max(yDEMCoordsTemp)]
                
        # deallocate arrays
        xDEMTemp = None
        yDEMTemp = None
        zDEMTemp = None
        xDEMResTemp = None
        yDEMResTemp = None
        xDEMCoordsTemp = None
        yDEMCoordsTemp = None
            
        # now find which edges are in which DEM's for processing later
        
        edgDict = {}
        
        totalEdgInfoTable = np.empty((numEdg,6))    # create empty array to hold dem values
        totalEdgInfoTable[:,0] = np.arange(numEdg)  # polulate edge numbers
        totalEdgInfoTable[:,1] = np.min(xS,axis=1)  # populate minimum x values of the vertices of the edge
        totalEdgInfoTable[:,2] = np.max(xS,axis=1)  # populate maximum x values of the vertices of the edge
        totalEdgInfoTable[:,3] = np.min(yS,axis=1)  # populate minimum y values of the vertices of the edge
        totalEdgInfoTable[:,4] = np.max(yS,axis=1)  # populate maximum y values of the vertices of the edge
        
        # find if each element is within any of the given polygons
        if hasattr(self.control, 'sgs_region_mask'):
            polys = gpd.read_file(self.control.sgs_region_mask).unary_union
            edgeInsidePolygon = \
            [all(gpd.GeoSeries([
                Point(xS[j,0], yS[j,0]),
                Point(xS[j,1], yS[j,1])]
                ).set_crs(epsg="4326").within(polys)) for j in range(numEdg)]
        else:
            edgeInsidePolygon = np.ones(numEdg).astype(bool)
            
        # loop through DEMs and determine which edges are in them
        # make sure to have DEMs in priority order meaning if you want to use fine
        # resolution in some areas have those dems listed first
        
        containedEdgList0Index = []
        noEdgDEMs = []
        
        for i in range(nDEMFilenameList):

            # find which edges are in the DEM
            totalEdgInfoTable[:,5] = ((totalEdgInfoTable[:,1]>=elevationDict["bounds%s"%i][0])
                                  & ((totalEdgInfoTable[:,2])<=elevationDict["bounds%s"%i][1])
                                  & ((totalEdgInfoTable[:,3])>=elevationDict["bounds%s"%i][2])
                                  & ((totalEdgInfoTable[:,4])<=elevationDict["bounds%s"%i][3]))
        
            whichAreInside = list(np.where([totalEdgInfoTable[j,5] == 1 and edgeInsidePolygon[j] for j in range(numEdg)])[0])
            edgDict["DEM%s"%i] = totalEdgInfoTable[whichAreInside,0].astype(int)        # store edge numbers of the edges inside the DEM bound
            
            whichAreInsideActualEdgNumber = totalEdgInfoTable[whichAreInside,0].astype(int) # get the actual edge numbers 
            totalEdgInfoTable = np.delete(totalEdgInfoTable,whichAreInside,axis=0)          # delete edges so we can skip those in the next dem

            # keep track if a dem does not have any elements inside to throw and exception later
            if len(whichAreInside) == 0:
                noEdgDEMs.append(self.control.demFilenameList[i])
            
            containedEdgList0Index.append(whichAreInsideActualEdgNumber)    # create a list of elements within subgrid area
        
        # throw exception if a dem has no element and print those dem names
        if(len(noEdgDEMs) != 0):
            for demName in noEdgDEMs:
                print(demName)
            sys.exit('No edges in the above DEMs, throw those puppies out\n and their matching landcover!\n')
                        
        # now concatenate the lists from above
        containedEdgList0Index = np.hstack(containedEdgList0Index)
        
        # now delete double counted vertices and edges
        containedEdgList0Index = np.unique(containedEdgList0Index).astype(int)
        
        # now I want to create a list of 1s and 0s to show whether or not a
        # vertex or edge is in the subgrid region
        
        binaryEdgList = np.zeros(numEdg)
        binaryEdgList[containedEdgList0Index] = 1
        
        # make an int array
        binaryEdgList = binaryEdgList.astype(int)
                
        # pre allocate arrays for subgrid quantities
        surfElevs = np.array([],dtype=np.float32)
        wetFraction = np.array([],dtype=np.float32)
        length = np.zeros((numEdg)).astype(np.float32)
        totWatDepth = np.array([],dtype=np.float32)
        wetTotWatDepth = np.array([],dtype=np.float32)
        cf = np.array([],dtype=np.float32)
        minElevationEdg = np.zeros(numEdg).astype(np.float32)    # lowest elevation in each edge for use in variable phi
        maxElevationEdg = np.zeros(numEdg).astype(np.float32)    # highest elevation in each edge for use in variable phi

        # fill length array
        length[np.where(binaryEdgList == 0)] = -99999
        
        # fill min/max elevation
        minElevationEdg[np.where(binaryEdgList == 0)] = 99999
        maxElevationEdg[np.where(binaryEdgList == 0)] = -99999

        # allocate the index that maps edge to location in vector
        edgIndex = np.zeros((numEdg)).astype(int)
        edgNumLevel = np.zeros((numEdg)).astype(int)
        edgIndexCnt = 0

        # create variable to keep track of what DEM you have read in
        
        # first create a loop for DEMs
        for i in range(nDEMFilenameList):
            
            # reading in DEM again
            # all variables the same as before
            if self.control.demFilenameList:
                elevationData = scm.importDEM(self.control.demFilenameList[i])
                xDEM = elevationData[0]
                yDEM = elevationData[1]
                zDEM = elevationData[2]
                xDEMRes = elevationData[3]
                yDEMRes = -1*elevationData[4]
                xDEMCoords = elevationData[5]
                yDEMCoords = elevationData[6]
                elevationData = None # deallocate 
            else:
                if not 'elevationData' in locals():
                    elevationData = self.importDEM_external()

                xDEM = elevationData[0]
                yDEM = elevationData[1]
                zDEM = elevationData[2]
                xDEMRes = elevationData[3]
                yDEMRes = -1*elevationData[4]
                xDEMCoords = elevationData[5]
                yDEMCoords = elevationData[6]
                elevationData = None # deallocate 

            # get a list of edges within this DEM
            edgList = edgDict["DEM%s"%i].astype(int)
            
            countEdgLoop = 0

            # loop through the edges
            for edg in edgList:
                startTime = time.perf_counter()
            
                # get x and y points for vertices of this element
                currXPerimeterPoints = xS[edg,:]
                currYPerimeterPoints = yS[edg,:]
                
                countInEdg = 0
                
                # get the bounds of the edge
                xyArray = np.array((currXPerimeterPoints,currYPerimeterPoints)).T

                # convert to meters for area calculations
                if self.need_projection:
                    xyCurrElementMeters = scm.projectMeshToMercator(xyArray[:,1],
                                                                    xyArray[:,0])
                else:
                    xyCurrElementMeters = [xyArray[:,0],xyArray[:,1]]
                
                # calculate the length of the edge
                edgLen = scm.seglength(xyCurrElementMeters[0][0],
                                       xyCurrElementMeters[1][0],
                                       xyCurrElementMeters[0][1],
                                       xyCurrElementMeters[1][1])
    
                # store edge length to array
                length[edg] = edgLen
            
                # cut down DEM further to save in computation expense
                
                # max and min X and Y for the edge
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
            
                # mask to find which cells are along an element edge
                mask = scm.isAlong(currXPerimeterPoints[0],currYPerimeterPoints[0],
                                   currXPerimeterPoints[1],currYPerimeterPoints[1],
                                   xCutGeoTiffMatrix2,yCutGeoTiffMatrix2,
                                   xDEMRes,yDEMRes)
                                
                # convert mask to cupy array
                mask = cp.asarray(mask)
                
                # count how many cells are within element
                countIn = cp.count_nonzero(mask)
                
                # if there are no cells within the element the DEM is too coarse
                # you must decrease the DEM resolution in this area
                if countIn == 0:
                    sys.exit('No DEM point along edge {}!'.format(edg))

                # keep track of this for use later
                countInEdg += countIn            
                
                zCutGeoTiffMatrix2 = cp.asarray(zCutGeoTiffMatrix2)
                zCutGeoTiffMatrix2masked = zCutGeoTiffMatrix2[mask]
                
                # get the min/max elevation of the element
                minElev = cp.nanmin(zCutGeoTiffMatrix2masked)
                maxElev = cp.nanmax(zCutGeoTiffMatrix2masked)
            
                # create array of surface elevations
                # this is used to calculate the subgrid variables for varying water depths
                surfElevIncrement = self.surfElevIncrement
                minElevFloor = math.floor(minElev/surfElevIncrement)*surfElevIncrement
                surfaceElevations = np.arange(minElevFloor,maxElev+surfElevIncrement,surfElevIncrement)
                if len(surfaceElevations) < 2:
                    surfaceElevations = np.array([minElev,maxElev]).astype(np.float32)
                else:
                    surfaceElevations[0] = minElev
                    surfaceElevations[-1] = maxElev
                r1 = maxElev + self.watDepthAboveHighestGroundIncrement
                r2 = r1 + self.maxWatDepthAboveHighestGround
                surfaceElevations = np.append(surfaceElevations, np.arange(r1,r2))
                num_SfcElevs = len(surfaceElevations)   # number of surface elevations we are running
                surfaceElevations = cp.asarray(surfaceElevations)

                # preallocate arrays
                wetDryEdgList = cp.zeros(num_SfcElevs)   # for wet area fraction phi
                totWatEdgList = cp.zeros(num_SfcElevs)   # for grid total water depth H_G
                
                ################ BEGIN VECTORIZED CALCULATIONS ##############################
                # create a 3d surface array for use in calculations
                tempSurfaceElevArray = cp.ones((len(zCutGeoTiffMatrix2masked),
                                               num_SfcElevs))*surfaceElevations
            
                # subtract the bathymetry (2D Array) array from the surface 
                # elevations (3D array) to get total water depths along the 
                # edge
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
                wetDryEdgList += tempCountWet

                ###############################################################

                #################### CALCULATING TOTAL WATER DEPTH ############
            
                # 0 out any dry cells
                tempTotWatDepthWetArray = tempTotWatDepthArray * tempWetDryList
                
                # integrate the wet total water depths for use in averaging later
                totWatEdgList += cp.sum(tempTotWatDepthWetArray,axis=0)

                ###############################################################
                            
                # Okay now we can finalize the values
                wetAvgTotWatDepth = totWatEdgList/wetDryEdgList
                gridAvgTotWatDepth = totWatEdgList/countInEdg
                wetFractionTemp = wetDryEdgList/countInEdg

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

                ####### store the values #######
                # store the index
                edgIndex[edg] = edgIndexCnt
                edgNumLevel[edg] = len(surfaceElevations)
                edgIndexCnt = edgIndexCnt + len(surfaceElevations)

                # store values
                surfElevs = np.append(surfElevs,surfaceElevations)
                totWatDepth = np.append(totWatDepth,gridAvgTotWatDepth)
                wetTotWatDepth = np.append(wetTotWatDepth,wetAvgTotWatDepth)
                wetFraction = np.append(wetFraction,wetFractionTemp)

                minElevationEdg[edg] = minElev
                maxElevationEdg[edg] = maxElev
                
                countEdgLoop += 1
                if countEdgLoop%1000==0:
                    stopTime = time.perf_counter()
                    print("Finished Edge {0} of {1} in DEM {2} took {3}".format(countEdgLoop,len(edgList),i,stopTime - startTime))
            
        wetFraction[np.isnan(wetFraction)] = 0.0
        wetTotWatDepth[np.isnan(wetTotWatDepth)] = 0.0
        totWatDepth[np.isnan(totWatDepth)] = 0.0

        # sort the vectors in the ascending order of the edg no.
        print('Sorting...',end='')
        def sortarray(vals):
            sorted = np.zeros(len(vals)).astype(np.float32)
            sortedIndexStart = np.zeros(numEdg).astype(int)
            sortedIndexEnd = np.zeros(numEdg).astype(int)
            sortedIndexCnt = 0
            for edg in range(numEdg):
                numLevels = edgNumLevel[edg]
                s = edgIndex[edg]
                e = edgIndex[edg]+numLevels
                snew = sortedIndexCnt
                enew = sortedIndexCnt+numLevels
                sorted[snew:enew] = vals[s:e]
                sortedIndexStart[edg] = snew
                sortedIndexEnd[edg] = enew
                sortedIndexCnt = enew
            return sorted, sortedIndexStart

        surfElevs, edgIndexSorted = sortarray(surfElevs)
        wetFraction, edgIndexSorted = sortarray(wetFraction)
        totWatDepth, edgIndexSorted = sortarray(totWatDepth)

        print('Done.')

        # Find global min/max surface elevations
        minElevationGlobal = np.min(minElevationEdg)
        maxElevationGlobal = np.max(maxElevationEdg)

        # store the resulting values
        self.subgridvectorized.edgIndex = edgIndexSorted
        self.subgridvectorized.surfaceElevationsEdg = surfElevs
        self.subgridvectorized.wetFractionEdg = wetFraction
        self.subgridvectorized.edglength = length
        self.subgridvectorized.totWatDepthEdg = totWatDepth
        self.subgridvectorized.binaryElementListEdg = binaryEdgList
        self.subgridvectorized.minElevationEdg = minElevationEdg
        self.subgridvectorized.maxElevationEdg = maxElevationEdg
        self.subgridvectorized.minElevationGlobal = minElevationGlobal
        self.subgridvectorized.maxElevationGlobal = maxElevationGlobal
        self.subgridvectorized.loadedEdg = True
        
    ############## WRITE SUBGRID CORRECTION DATA TO NETCDF FOR VECTORIZED STORAGE ##############

    def writeSubgridLookupTableNetCDFForVectorizedStorage(self):
        from netCDF4 import Dataset
        import numpy as np

        # open netcdf file
        ncFile = Dataset(self.control.outputFilename, mode = 'w', format = 'NETCDF4')
        
        # Set version number
        ncFile.data_format = 'v1.0.0'

        # create dimensions
        ncFile.createDimension('numElem',self.mesh.numEle)                                    # element dimension
        ncFile.createDimension('numNode',self.mesh.numNod)                                    # number of nodes in mesh
        ncFile.createDimension('numEdge',self.mesh.numEdg)                                    # edge dimension
        ncFile.createDimension('numEdgePerElem',3)                                            # edge dimension
        ncFile.createDimension('vecLen',len(self.subgridvectorized.surfaceElevations))    # vector length
        ncFile.createDimension('edgeVecLen',len(self.subgridvectorized.surfaceElevationsEdg)) # vector length
        
        # create variables
        # write variable dimensions transposed because FORTRAN will read them that way
        # only need to do it for the 3D arrays, handle 2 and 1D a different way.
        # - element variables
        elemIndex = ncFile.createVariable('elemLocations',int,'numElem')                       # indexes showing the locations of elements in the vectorized storage        
        wetFractionVarElement = ncFile.createVariable('wetFraction',np.float32,'vecLen')  # elemental wet area fraction
        areaVar = ncFile.createVariable('area',np.float32,'numElem')                              # elemental areas
        totWatDepthVar = ncFile.createVariable('totWatDepth',np.float32,'vecLen')         # elemental grid averaged total water depth
        surfaceElevationsVar = ncFile.createVariable('surfaceElevations',np.float32,'vecLen') # surface elevation array
        cfVarElement = ncFile.createVariable('cf',np.float32,'vecLen')                    # elemental coefficient of friction level 0
        if self.level0andLevel1:                                                                  # elemental coefficient of friction level 1
            cmfVarElement = ncFile.createVariable('cmf',np.float32,'vecLen')
        binaryElementListVariable = ncFile.createVariable('binaryList',int,'numElem')      # variables showing which elements are contained within the subgrid area
        minElevationEleVariable = ncFile.createVariable('minElevation',np.float32,'numElem')  # min elevation
        maxElevationEleVariable = ncFile.createVariable('maxElevation',np.float32,'numElem')  # max elevation
        # - edge variables
        edgeIndex = ncFile.createVariable('edgeLocations',int,'numEdge')                       # indexes showing the locations of elements in the vectorized storage        
        lenVarEdge = ncFile.createVariable('edgeLength',np.float32,'numEdge')                     # elemental areas
        totWatDepthVarEdge = ncFile.createVariable('edgeTotWatDepth',np.float32,'edgeVecLen')     # elemental grid averaged total water depth
        surfaceElevationsVarEdge = ncFile.createVariable('edgeSurfaceElevations',np.float32,'edgeVecLen') # surface elevation array
        binaryEdgeListVariable = ncFile.createVariable('edgeBinaryList',int,'numEdge')         # variables showing which elements are contained within the subgrid area
        minElevationEdgeVariable = ncFile.createVariable('edgeMinElevation',np.float32,'numEdge') # min elevation
        maxElevationEdgeVariable = ncFile.createVariable('edgeMaxElevation',np.float32,'numEdge') # max elevation
        # - element-to-edge table
        el2ed = ncFile.createVariable('el2ed',int,('numElem','numEdgePerElem'))                # table to relate an element to its edges        
        
        # assgin values
        # - element values
        elemIndex[:] = self.subgridvectorized.elemIndex + 1                                       # shift elemIndex to make it begin with 1
        wetFractionVarElement[:] = self.subgridvectorized.wetFraction
        areaVar[:] = self.subgridvectorized.area
        totWatDepthVar[:] = self.subgridvectorized.totWatDepth
        surfaceElevationsVar[:] = self.subgridvectorized.surfaceElevations
        cfVarElement[:] = self.subgridvectorized.cf
        if self.level0andLevel1:
            cmfVarElement[:] = self.subgridvectorized.cmf
        binaryElementListVariable[:] = self.subgridvectorized.binaryElementList
        minElevationEleVariable[:] = self.subgridvectorized.minElevationEle
        maxElevationEleVariable[:] = self.subgridvectorized.maxElevationEle
        # - edge values
        edgeIndex[:] = self.subgridvectorized.edgIndex + 1                                       # shift edgeIndex to make it begin with 1
        lenVarEdge[:] = self.subgridvectorized.edglength
        totWatDepthVarEdge[:] = self.subgridvectorized.totWatDepthEdg
        surfaceElevationsVarEdge[:] = self.subgridvectorized.surfaceElevationsEdg
        binaryEdgeListVariable[:] = self.subgridvectorized.binaryElementListEdg
        minElevationEdgeVariable[:] = self.subgridvectorized.minElevationEdg
        maxElevationEdgeVariable[:] = self.subgridvectorized.maxElevationEdg
        # - element-to-edge table
        el2ed[:,:] = self.mesh.el2ed + 1

        # close netcdf file
        ncFile.close()


    ######## READ SUBGRID CORRECTION DATA TO NETCDF FOR VECTORIZED STORAGE ########

    def readSubgridLookupTableNetCDFForVectorizedStorage(self):
        import numpy as np
        from netCDF4 import Dataset

        # open netcdf file
        lookupTable = Dataset(self.control.outputFilename)

        # read variables
        # - element variables
        self.subgridvectorized.elemIndex = np.asarray(lookupTable['elemLocations'][:])-1
        self.subgridvectorized.surfaceElevations = np.asarray(lookupTable['surfaceElevations'][:])
        self.subgridvectorized.wetFraction = np.asarray(lookupTable['wetFraction'][:])
        self.subgridvectorized.area = np.asarray(lookupTable['area'][:])
        self.subgridvectorized.totWatDepth = np.asarray(lookupTable['totWatDepth'][:])
        self.subgridvectorized.cf = np.asarray(lookupTable['cf'][:]).T
        if self.level0andLevel1:
            self.subgridvectorized.cmf = np.asarray(lookupTable['cmf'][:]).T
        self.subgridvectorized.binaryElementList = np.asarray(lookupTable['binaryList'][:])
        self.subgridvectorized.minElevationEle = np.asarray(lookupTable['minElevation'][:])
        self.subgridvectorized.maxElevationEle = np.asarray(lookupTable['maxElevation'][:])
        # - edge variables
        self.subgridvectorized.edgIndex = np.asarray(lookupTable['edgeLocations'][:])-1
        self.subgridvectorized.surfaceElevationsEdg = np.asarray(lookupTable['edgeSurfaceElevations'][:])
        self.subgridvectorized.edglength = np.asarray(lookupTable['edgeLength'][:])
        self.subgridvectorized.totWatDepthEdg = np.asarray(lookupTable['edgeTotWatDepth'][:])
        self.subgridvectorized.binaryEdgList = np.asarray(lookupTable['edgeBinaryList'][:])
        self.subgridvectorized.minElevationEdg = np.asarray(lookupTable['edgeMinElevation'][:])
        self.subgridvectorized.maxElevationEdg = np.asarray(lookupTable['edgeMaxElevation'][:])
        # - element-to-edge table
        self.subgridvectorized.el2ed = np.asarray(lookupTable['el2ed'][:,:])-1

        lookupTable.close()

    ######## Plot subgrid correction data for vectorized storage ########

    def plot_subgrid_for_vectorizedstorage(self,imagePath):
        import numpy as np
        from matplotlib.patches import Circle, Wedge, Polygon
        from matplotlib.collections import PatchCollection
        import matplotlib.pyplot as plt
        import os
 
        # Check whether the image path exists or not and create it if not
        isExist = os.path.exists(imagePath)
        if not isExist:
            os.makedirs(imagePath)
   
        # create polygons
        patches = []
        for i in range(self.mesh.numEle):
            t = self.mesh.tri[i,:]
            x = np.asarray(self.mesh.coord['Longitude'])[t]
            y = np.asarray(self.mesh.coord['Latitude'])[t]
            p = np.concatenate((x[:,np.newaxis],y[:,np.newaxis]), axis=1)

            polygon = Polygon(p, closed=False)
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
            ax.set_aspect('equal','box')
            # plt.show()
            fig.savefig("{:s}/sg_{:s}.png".format(imagePath,figfile),dpi=600,facecolor='white',edgecolor='none')
            plt.close()

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
                    if se < minElevationEle[ele]:
                        vinterp[ele] = np.nan
                    elif  se > maxElevationEle[ele]:
                        if ele == numEle - 1:
                            ii = len(ve) - 1
                        else:
                            ii = elemIndex[ele+1] - 1
                        vinterp[ele] = ve[ii]
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
                            if ista == iend:
                                break
                        if found == False:
                            print("SE:{:.2f}, surfaceElev[0,end]=({:.2f},{:.2f})".format(se,surfaceEle[0],surfaceEle[-1]))
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
                ax.set_aspect('equal','box')
                # plt.show()
                fig.savefig("{:s}/sg_{:s}_{:03d}.png".format(imagePath,figfile,i),dpi=600,facecolor='white',edgecolor='none')
                plt.close()

        surfElevForPlot = np.arange(-20.0,0.0,4.0)

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


