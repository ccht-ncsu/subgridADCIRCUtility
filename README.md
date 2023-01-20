# subgridADCIRCUtility
This set of scripts will create subgrid input files for subgrid enabled ADCIRC.

## Example:

The following link will take you to a data repository (need to figure out where to put this) where you will find an example folder containing a text control file and a python script used to call the subgrid calculator for a ADCIRC mesh of Galveston Bay, TX along with several other datafiles needed to perform the subgrid calulations.

1. The control file lists:
  - The outputfile name for the subgrid lookup table.
  - The ADCIRC mesh file.
  - A flag saying whether or not you are using default mannings n values.
  - If not using default mannings n values you need to add a line with the filename of the non-defualt table.
  - The number of DEMs used in the calulation.
  - The DEM file(s) used in the subgrid calulation.
  - The number of landcover files used in the calulation.
  - The landcover file(s) used in the calulation.

2. The python script reads the control file, calls the subgrid calulator and then simplifies the lookup table.
  - When calling the subgrid calculator there are 2 versions to choose from. One is a CPU only version, and the other uses both CPU and NVIDIA GPUs. The CPU/GPU code         will run faster, especially for high-resolution subgrid files and large meshes.
  - The GPU version of the code is called from the function calculateSubgridCorrectionv3GPU, and the CPU version is called using calculateSubgridCorrectionv3.

  
3. Other files contained within the data repo include:

  - 10 m resolution DEM and Landcover .tif files of Galevston Bay.
  - A relatively coarse ADCIRC mesh file, nodal attribute file, and fort.15 file for use with subgrid ADCIRC.
  - Some visualizations of subgrid ADCIRC and Conventional ADCIRC run on this domain.
  - Fully completed subgrid NetCDF lookup files to compare your results to.
  - A mannings n table with non-default values.
  
## PYTHON MODULES NEEDED

1. pandas
2. numpy
3. matplotlib
4. gdal
5. time
6. netCDF4
7. cupy (for GPU)
8. cmocean (if you want to visualize stuff)

# Gallery

This is an example of the vertex averaged variables produced for Galveston Bay, Texas.
![GBAY_git_figure_fixed_again](https://user-images.githubusercontent.com/50885561/189756706-9c3d6b2c-af4b-4725-a3b6-009e42d8e340.jpeg)


