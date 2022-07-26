# subgridADCIRCUtility
This set of scripts will create subgrid input files for subgrid enabled ADCIRC.

Example:

The following link will take you to a data repository (need to figure out where to put this) where you will find an example folder containing a text control file and a python script used to call the subgrid calculator for a ADCIRC mesh of Galveston Bay, TX along with several other datafiles needed to perform the subgrid calulations.

1. The control file lists the outputfile name for the subgrid lookup table, the ADCIRC mesh file, the DEM file used in the subgrid calulations, and the landcover file also used in the calulations.

2. The python script reads the control file, calls the subgrid calulator and then simplifies the lookup table.
  - When calling the subgrid calculator there are 2 versions to choose from. One is a CPU only version, and the other uses both CPU and NVIDIA GPUs. The CPU/GPU code         will run faster, especially for high-resolution subgrid files and large meshes. 
  - There is a separate function call to reduce the lookup table. You HAVE to do this if you want to run subgrid ADCIRC. I have not incorporated this reduction function     in the main calulator script yet but will soon.
  
3. Other files contained within the data repo include:

  - 10 m resolution DEM and Landcover .tif files of Galevston Bay.
  - A relatively coarse ADCIRC mesh file, nodal attribute file, and fort.15 file for use with subgrid ADCIRC.
  - Some visualizations of subgrid ADCIRC and Conventional ADCIRC run on this domain.
  - Fully completed subgrid NetCDF lookup files to compare your results to.
  
PYTHON MODULES NEEDED

1. pandas
2. numpy
3. matplotlib
4. gdal
5. time
6. netCDF4
7. cupy (for GPU)
8. cmocean (if you want to visualize stuff)
