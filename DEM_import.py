import numpy as np

def importDEM_QuarterAnnularPolynomial():
    print('Pseudo DEM generation process started.')

    ft2m = 0.3048
    buffDist = 10000.0
    xMin = -buffDist
    xMax = 152400.0 + buffDist
    yMin = xMin
    yMax = xMax
    xDim = 3001
    yDim = 3001
    xRes = + (xMax - xMin) / float(xDim-1)
    yRes = - (yMax - yMin) / float(yDim-1)
    x_upper_left = xMin
    y_upper_left = yMax
    
    xCoords = np.zeros((xDim))
    yCoords = np.zeros((yDim))
    
    for i in range(xDim):
        xCoords[i] = x_upper_left + (xRes)*i
    
    for i in range(yDim):
        yCoords[i] = y_upper_left + (yRes)*i
    
    [X,Y] = np.meshgrid(xCoords,yCoords)

    z_array = np.zeros((yDim,xDim))

    for j in range(yDim):
        if j%1000 == 0:
            print(" - now at j = {} / {}".format(j,yDim))
        for i in range(xDim):
            x = xCoords[i]
            y = yCoords[j]
            r = np.sqrt(x*x + y*y) / ft2m
            z_array[j,i] = - 2.5e-10 * r * r * ft2m

    print('Pseudo DEM generation process finished.')

    return X,Y,z_array,xRes,yRes,xCoords,yCoords

def importDEM_QuarterAnnularPolynomialSmallDomain():
    print('Pseudo DEM generation process started.')

    hscale = 0.01
    ft2m = 0.3048
    buffDist = 10000.0 * hscale
    xMin = -buffDist
    xMax = 152400.0 * hscale + buffDist
    yMin = xMin
    yMax = xMax
    xDim = 3001
    yDim = 3001
    xRes = + (xMax - xMin) / float(xDim-1)
    yRes = - (yMax - yMin) / float(yDim-1)
    x_upper_left = xMin
    y_upper_left = yMax
    
    xCoords = np.zeros((xDim))
    yCoords = np.zeros((yDim))
    
    for i in range(xDim):
        xCoords[i] = x_upper_left + (xRes)*i
    
    for i in range(yDim):
        yCoords[i] = y_upper_left + (yRes)*i
    
    [X,Y] = np.meshgrid(xCoords,yCoords)

    z_array = np.zeros((yDim,xDim))

    for j in range(yDim):
        if j%1000 == 0:
            print(" - now at j = {} / {}".format(j,yDim))
        for i in range(xDim):
            x = xCoords[i]
            y = yCoords[j]
            r = np.sqrt(x*x + y*y) / ft2m / hscale
            z_array[j,i] = - 2.5e-10 * r * r * ft2m

    print('Pseudo DEM generation process finished.')

    return X,Y,z_array,xRes,yRes,xCoords,yCoords

def importDEM_QuarterAnnularPolynomialSmallDomainWetDry():
    print('Pseudo DEM generation process started.')

    hscale = 0.01
    zshift = 4.0
    ft2m = 0.3048
    buffDist = 10000.0 * hscale
    xMin = -buffDist
    xMax = 152400.0 * hscale + buffDist
    yMin = xMin
    yMax = xMax
    xDim = 3001
    yDim = 3001
    xRes = + (xMax - xMin) / float(xDim-1)
    yRes = - (yMax - yMin) / float(yDim-1)
    x_upper_left = xMin
    y_upper_left = yMax
    
    xCoords = np.zeros((xDim))
    yCoords = np.zeros((yDim))
    
    for i in range(xDim):
        xCoords[i] = x_upper_left + (xRes)*i
    
    for i in range(yDim):
        yCoords[i] = y_upper_left + (yRes)*i
    
    [X,Y] = np.meshgrid(xCoords,yCoords)

    z_array = np.zeros((yDim,xDim))

    for j in range(yDim):
        if j%1000 == 0:
            print(" - now at j = {} / {}".format(j,yDim))
        for i in range(xDim):
            x = xCoords[i]
            y = yCoords[j]
            r = np.sqrt(x*x + y*y) / ft2m / hscale
            z_array[j,i] = - 2.5e-10 * r * r * ft2m + zshift

    print('Pseudo DEM generation process finished.')

    return X,Y,z_array,xRes,yRes,xCoords,yCoords

def importDEM_QuarterAnnularPolynomialSmallDomainWetDry2():
    print('Pseudo DEM generation process started.')

    hscale = 0.01
    zshift = 10.0
    ft2m = 0.3048
    buffDist = 10000.0 * hscale
    xMin = -buffDist
    xMax = 152400.0 * hscale + buffDist
    yMin = xMin
    yMax = xMax
    xDim = 3001
    yDim = 3001
    xRes = + (xMax - xMin) / float(xDim-1)
    yRes = - (yMax - yMin) / float(yDim-1)
    x_upper_left = xMin
    y_upper_left = yMax
    
    xCoords = np.zeros((xDim))
    yCoords = np.zeros((yDim))
    
    for i in range(xDim):
        xCoords[i] = x_upper_left + (xRes)*i
    
    for i in range(yDim):
        yCoords[i] = y_upper_left + (yRes)*i
    
    [X,Y] = np.meshgrid(xCoords,yCoords)

    z_array = np.zeros((yDim,xDim))

    for j in range(yDim):
        if j%1000 == 0:
            print(" - now at j = {} / {}".format(j,yDim))
        for i in range(xDim):
            x = xCoords[i]
            y = yCoords[j]
            r = np.sqrt(x*x + y*y) / ft2m / hscale
            z_array[j,i] = - 2.5e-10 * r * r * ft2m + zshift

    print('Pseudo DEM generation process finished.')

    return X,Y,z_array,xRes,yRes,xCoords,yCoords

def importDEM_QuarterAnnularPolynomialRectSmallDomain():
    print('Pseudo DEM generation process started.')

    hscale = 0.01
    ft2m = 0.3048
    buffDist = 10000.0 * hscale
    xMin = -buffDist
    xMax = 152400.0 * hscale + buffDist
    yMin = -buffDist
    yMax = 152400.0 * hscale * 4.0 / 3.0 + buffDist
    xDim = 3001
    yDim = int(3001 * 4.0 / 3.0)
    xRes = + (xMax - xMin) / float(xDim-1)
    yRes = - (yMax - yMin) / float(yDim-1)
    x_upper_left = xMin
    y_upper_left = yMax
    
    xCoords = np.zeros((xDim))
    yCoords = np.zeros((yDim))
    
    for i in range(xDim):
        xCoords[i] = x_upper_left + (xRes)*i
    
    for i in range(yDim):
        yCoords[i] = y_upper_left + (yRes)*i
    
    [X,Y] = np.meshgrid(xCoords,yCoords)

    z_array = np.zeros((yDim,xDim))

    for j in range(yDim):
        if j%1000 == 0:
            print(" - now at j = {} / {}".format(j,yDim))
        for i in range(xDim):
            x = xCoords[i]
            r = x / ft2m / hscale
            z_array[j,i] = - 2.5e-10 * r * r * ft2m

    print('Pseudo DEM generation process finished.')

    return X,Y,z_array,xRes,yRes,xCoords,yCoords

def importDEM_QuarterAnnularPolynomialRectSmallDomainFlat():
    print('Pseudo DEM generation process started.')

    hscale = 0.01
    ft2m = 0.3048
    buffDist = 10000.0 * hscale
    xMin = -buffDist
    xMax = 152400.0 * hscale + buffDist
    yMin = -buffDist
    yMax = 152400.0 * hscale * 4.0 / 3.0 + buffDist
    xDim = 3001
    yDim = int(3001 * 4.0 / 3.0)
    xRes = + (xMax - xMin) / float(xDim-1)
    yRes = - (yMax - yMin) / float(yDim-1)
    x_upper_left = xMin
    y_upper_left = yMax
    
    xCoords = np.zeros((xDim))
    yCoords = np.zeros((yDim))
    
    for i in range(xDim):
        xCoords[i] = x_upper_left + (xRes)*i
    
    for i in range(yDim):
        yCoords[i] = y_upper_left + (yRes)*i
    
    [X,Y] = np.meshgrid(xCoords,yCoords)

    z_array = np.zeros((yDim,xDim))

    for j in range(yDim):
        if j%1000 == 0:
            print(" - now at j = {} / {}".format(j,yDim))
        for i in range(xDim):
            x = xCoords[i]
            r = x / ft2m / hscale
            z_array[j,i] = - 19.0

    print('Pseudo DEM generation process finished.')

    return X,Y,z_array,xRes,yRes,xCoords,yCoords

