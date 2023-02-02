# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 18:36:30 2023

@author: jlwoodr3
"""

# this script will play with reducing the vertex lookup tables

import sys 
import numpy as np
import matplotlib.pyplot as plt
import cmocean
import netCDF4 as nc
import matplotlib as mpl
import matplotlib.ticker as mtick

sys.path.append(r'C:\Users\jlwoodr3\Documents\GitHub\subgridADCIRCUtility')

import subgrid_calculator as sc
# first read in the lookup table to practice on

lookupTableFilename = r'C:\Users\jlwoodr3\Desktop\Oceanmesh_data_and_scripts\subgrid_GBAY_practice_vectorize\subgrid\avgVars_02_increment_dev.nc'

lookupTable = nc.Dataset(lookupTableFilename)

# load in all of the vertex variables and the phi set

phiSet = np.asarray(lookupTable['phiSet'][:])

wetFractionVertex = np.asarray(lookupTable['wetFractionVertex'][:])

gridTotWatDepthVertex = np.asarray(lookupTable['gridTotWatDepthVertex'][:])

wetTotWatDepthVertex = np.asarray(lookupTable['wetTotWatDepthVertex'][:])

cf = np.asarray(lookupTable['cfVertex'][:])

cmf = np.asarray(lookupTable['cmfVertex'][:])

surfaceElevations = np.asarray(lookupTable['surfaceElevations'][:])

# now practice redution

depthsVertForLookup = np.zeros((11,len(wetFractionVertex[:])))
HGVertForLookup = np.zeros((11,len(wetFractionVertex[:])))
HWVertForLookup = np.zeros((11,len(wetFractionVertex[:])))
cfVertForLookup = np.zeros((11,len(wetFractionVertex[:])))
cmfVertForLookup = np.zeros((11,len(wetFractionVertex[:])))

# # first fill the elements that are either all dry or all wet and we will
# # interpolate the ones in between
# # go ahead and set any subelement that does not get fully wet to dry
# checkWhere = np.all(wetFractionVertex<1.0,axis=1).nonzero()
# depthsVertForLookup[:,checkWhere] = np.max(surfaceElevations)
# HGVertForLookup[:,checkWhere] = 0.0
# HWVertForLookup[:,checkWhere] = 0.0
# cfVertForLookup[:,checkWhere] = np.max(cf[checkWhere,:])
# cmfVertForLookup[:,checkWhere] = np.max(cmf[checkWhere,:])
# # go ahead and set any element that is always wet to wet and we will interpolate later
# checkWhere = np.all(wetFractionVertex==1.0,axis=1).nonzero()
# depthsVertForLookup[:,checkWhere] = np.min(surfaceElevations)
# HGVertForLookup[:,checkWhere] = np.min(gridTotWatDepthVertex[checkWhere,:])
# HWVertForLookup[:,checkWhere] = np.min(wetTotWatDepthVertex[checkWhere,:])
# cfVertForLookup[:,checkWhere] = np.min(cf[checkWhere,:])
# cmfVertForLookup[:,checkWhere] = np.min(cmf[checkWhere,:])
# # now find the partially wet elements 
# checkWhere = np.any(wetFractionVertex == 0.0,axis=1).nonzero()[0]

for i in range(len(wetFractionVertex[:,0])):
    
    vert = i
    currPhiArray = wetFractionVertex[vert,:]
    
    # make sure that the phi array also gets fully wet and then proceed
    # otherwise just skip
    
    # for phi == 0 you want to find exactly where that is in the currPhiArray
    equalTo0 = np.where(currPhiArray == 0.0)[0]
    
    if(len(equalTo0)!=0): # if 0.0 exists in the array
    
        depthsVertForLookup[0,vert] = surfaceElevations[equalTo0[-1]]
        HGVertForLookup[0,vert] = gridTotWatDepthVertex[vert,equalTo0[-1]]
        HWVertForLookup[0,vert] = wetTotWatDepthVertex[vert,equalTo0[-1]]
        cfVertForLookup[0,vert] = cf[vert,equalTo0[-1]]
        cmfVertForLookup[0,vert] = cmf[vert,equalTo0[-1]]
        
    else: # so if it never gets fully dry set everything to the value corresponding to the first surface elevations
    
        depthsVertForLookup[0,vert] = surfaceElevations[0]
        HGVertForLookup[0,vert] = gridTotWatDepthVertex[vert,0]
        HWVertForLookup[0,vert] = wetTotWatDepthVertex[vert,0]
        cfVertForLookup[0,vert] = cf[vert,0]
        cmfVertForLookup[0,vert] = cmf[vert,0]
        
    # now check for when phi == 1.0 and find exactly where that is
    
    equalTo1 = np.where(currPhiArray == 1.0)[0]
    
    if(len(equalTo1)!=0): # if 1.0 exists in the array
    
        depthsVertForLookup[-1,vert] = surfaceElevations[equalTo1[0]]
        HGVertForLookup[-1,vert] = gridTotWatDepthVertex[vert,equalTo1[0]]
        HWVertForLookup[-1,vert] = wetTotWatDepthVertex[vert,equalTo1[0]]
        cfVertForLookup[-1,vert] = cf[vert,equalTo1[0]]
        cmfVertForLookup[-1,vert] = cmf[vert,equalTo1[0]]
        
    else: # if there is nothing that is equal to 1 (so never gets fully wet, just set everything to correspind to the last surface elevation)
    
        depthsVertForLookup[-1,vert] = surfaceElevations[-1]
        HGVertForLookup[-1,vert] = gridTotWatDepthVertex[vert,-1]
        HWVertForLookup[-1,vert] = wetTotWatDepthVertex[vert,-1]
        cfVertForLookup[-1,vert] = cf[vert,-1]
        cmfVertForLookup[-1,vert] = cmf[vert,-1]
        
        
    # now for everything else
    
    for k in range(1,len(phiSet)-1):
    
        desiredPhi = phiSet[k]
        greaterThan = np.where(currPhiArray > desiredPhi)[0]
        
        if(len(greaterThan)==0): # so if the first currphi index is greater than the desired phi 
        
        # set everything to correspond to the first surfaceelevation
            
            depthsVertForLookup[k,vert] = surfaceElevations[-1]
            HGVertForLookup[k,vert] = gridTotWatDepthVertex[vert,-1]
            HWVertForLookup[k,vert] = wetTotWatDepthVertex[vert,-1]
            cfVertForLookup[k,vert] = cf[vert,-1]
            cmfVertForLookup[k,vert] = cmf[vert,-1]
            
        elif(greaterThan[0] == 0): # so if nothing in the currPhiArray is greater than the desired phi
        
        # set everything to correspond to the last surface elevation
        
            depthsVertForLookup[k,vert] = surfaceElevations[0]
            HGVertForLookup[k,vert] = gridTotWatDepthVertex[vert,0]
            HWVertForLookup[k,vert] = wetTotWatDepthVertex[vert,0]
            cfVertForLookup[k,vert] = cf[vert,0]
            cmfVertForLookup[k,vert] = cmf[vert,0]

            
        else: # this is where we interpolate 
            
            greaterThan = greaterThan[0]
            lessThan = greaterThan - 1
    
            
            depthsVertForLookup[k,vert] = (((desiredPhi - currPhiArray[lessThan])
                                          /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                          *(surfaceElevations[greaterThan] - surfaceElevations[lessThan])
                                          + (surfaceElevations[lessThan]))
            
            HGVertForLookup[k,vert] = (((desiredPhi - currPhiArray[lessThan])
                                          /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                          *(gridTotWatDepthVertex[vert,greaterThan]
                                            - gridTotWatDepthVertex[vert,lessThan])
                                          + (gridTotWatDepthVertex[vert,lessThan]))
            
            HWVertForLookup[k,vert] = (((desiredPhi - currPhiArray[lessThan])
                                          /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                          *(wetTotWatDepthVertex[vert,greaterThan]
                                            - wetTotWatDepthVertex[vert,lessThan])
                                          + (wetTotWatDepthVertex[vert,lessThan]))
            cfVertForLookup[k,vert] = (((desiredPhi - currPhiArray[lessThan])
                                          /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                          *(cf[vert,greaterThan]
                                            - cf[vert,lessThan])
                                          + (cf[vert,lessThan]))
            cmfVertForLookup[k,vert] = (((desiredPhi - currPhiArray[lessThan])
                                          /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
                                          *(cmf[vert,greaterThan]
                                            - cmf[vert,lessThan])
                                          + (cmf[vert,lessThan]))
    
    
    # if(1.0 in currPhiArray):
        
    #     # do when phi == 0
    #     # for phi == 0.0 you want exactly where that is in the currPhiArray
    #     equalTo = np.where(currPhiArray == 0.0)[0][-1]
    #     depthsVertForLookup[0,vert] = surfaceElevations[equalTo]
    #     HGVertForLookup[0,vert] = gridTotWatDepthVertex[vert,equalTo]
    #     HWVertForLookup[0,vert] = wetTotWatDepthVertex[vert,equalTo]
    #     cfVertForLookup[0,vert] = cf[vert,equalTo]
    #     cmfVertForLookup[0,vert] = cmf[vert,equalTo]
        
    #     # do when phi == 1
    #     # for phi == 1.0 you want exactly where that is in the currPhiArray
    #     equalTo = np.where(currPhiArray == 1.0)[0][0]
    #     depthsVertForLookup[-1,vert] = surfaceElevations[equalTo]
    #     HGVertForLookup[-1,vert] = gridTotWatDepthVertex[vert,equalTo]
    #     HWVertForLookup[-1,vert] = wetTotWatDepthVertex[vert,equalTo]
    #     cfVertForLookup[-1,vert] = cf[vert,equalTo]
    #     cmfVertForLookup[-1,vert] = cmf[vert,equalTo]
        
        
    #     # only iterate between the second and next to last phi
    #     for k in range(1,len(phiSet)-1):
        
    #         desiredPhi = phiSet[k]
    #         greaterThan = np.where(currPhiArray > desiredPhi)[0][0]
    #         lessThan = greaterThan - 1
            
    #         # if(k == 1):
                
    #         #     print('check')
            
    #         depthsVertForLookup[k,vert] = (((desiredPhi - currPhiArray[lessThan])
    #                                       /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
    #                                       *(surfaceElevations[greaterThan] - surfaceElevations[lessThan])
    #                                       + (surfaceElevations[lessThan]))
            
    #         HGVertForLookup[k,vert] = (((desiredPhi - currPhiArray[lessThan])
    #                                       /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
    #                                       *(gridTotWatDepthVertex[vert,greaterThan]
    #                                         - gridTotWatDepthVertex[vert,lessThan])
    #                                       + (gridTotWatDepthVertex[vert,lessThan]))
            
    #         HWVertForLookup[k,vert] = (((desiredPhi - currPhiArray[lessThan])
    #                                       /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
    #                                       *(wetTotWatDepthVertex[vert,greaterThan]
    #                                         - wetTotWatDepthVertex[vert,lessThan])
    #                                       + (wetTotWatDepthVertex[vert,lessThan]))
    #         cfVertForLookup[k,vert] = (((desiredPhi - currPhiArray[lessThan])
    #                                       /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
    #                                       *(cf[vert,greaterThan]
    #                                         - cf[vert,lessThan])
    #                                       + (cf[vert,lessThan]))
    #         cmfVertForLookup[k,vert] = (((desiredPhi - currPhiArray[lessThan])
    #                                       /(currPhiArray[greaterThan] - currPhiArray[lessThan]))
    #                                       *(cmf[vert,greaterThan]
    #                                         - cmf[vert,lessThan])
    #                                       + (cmf[vert,lessThan]))
            
            

# checkwhere1 = np.any(wetFraction==1.0,axis=0).nonzero()
# checkwhereEle1 = checkwhere1[1]
# checkwhereVert1  = checkwhere1[0]
# depthsEleForLookup[:,checkwhereVert1,checkwhereEle1] = minSurfElev
# HEleForLookup[:,checkwhereVert1,checkwhereEle1] = totWatDepth[0,checkwhereVert1,checkwhereEle1]
# cadvForLookup[:,checkwhereVert1,checkwhereEle1] = cadv[-1,checkwhereVert1,checkwhereEle1]
# # now find where there are dry to partially wet to fully wet subelement
# checkwhere0 = np.any(wetFraction == 0.0,axis=0).nonzero()
# checkwhereVert0 = checkwhere0[0]
# checkwhereEle0 = checkwhere0[1]

