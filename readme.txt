This database has four folders:

-----------------------------------------------------------------------------------------------------

1.CodeN2fix
This folder contians the code for Terrestrial Ecosystem Model (TEM) to simulation N2 fixation. 
The codes are written in C++ and need to be compiled under Linux environment. 
To compile, go to the code folder, run: mpiicpc xtem423e1.cpp -o xtem423e1
The intel and impi modules need to be loaded before compiling. 

-----------------------------------------------------------------------------------------------------

2.runexample
This is an example folder of how to run TEM on a relatively small area. 
The data are used to simulate the biophysical process of some pan-Arctic grids given the focring data and parameter files. 
The universal parameters are under folder run/, while the spatially-explicit inputs and parameters are under part-0/ and part-1/. 

2.1 In part-0/ and part-1/, the files include:
This folder contains the input and parameter to run TEM. 
The entire dataset used to simulate global N2 fixation is not provided since it is too large. 
In stead, we give sample input folders: part-0/ and part-1/.

climate inputs: 
each row represent for input for one year, and has 22 columns. These columns are:
longitude, latitude, variable name, grid cell land area (km2), year , total of Jan-Dec, maximum of Jan-Dec, mean of Jan-Dec, minimum of Jan-Dec, Jan-Dec (12 columns), and region. 
These inputs include:
	clds.txt-*: cloudiness in %.
	prec.txt-*: precipitation in mm.
	tair.txt-*: temperature in degree °C.
	vap.txt-*: vapor pressure in hPa.
        Ndeponhx.txt-*: Atmospheric NHx deposition in gN/m-2
        Ndeponoy.txt-*: Atmospheric NOy deposition in gN/m-2

spatially-explicit inputs:
	ph.txt-*: pH input. The columns are lon, lat, variable name, grid cell land area (km2), pH and region.
	elev.txt-*: elevation (m) input. The columns are lon, lat, variable name, elevation and region.
	clfao.txt-*: soil texture input. The columns are lon, lat, variable nuame, grid cell land area (km2), present of sand, silt and clay, source, region.
	veg.txt-*: vegetaiton class input. The columns are lon, lat, variable name, dummy, vegetation class, region, region.
        density.txt-*: Soil bulk density (g cm-3).The columns are lon, lat, variable name, soil bulk density and region.

Please note that although some inputs have grid cell land area column, this value is not used within TEM. 
Therefore, even the values do not match in different files, it won't influence the simulation results. 

2.2 Under runexample folder, there are parameter files universal to all grid cells.
kco21860.data: co2 concentration (ppm).
run.go4: input file list of TEM. 
QTLA44A2.ECD: soil thermal properties in each layer.
QTSP44A2.ECD: soil thermal properties not specified for layers.
QTST44A2.ECD: initial soil thermal profile. 
Tcomm423.ecd: information of vegetation mosaic.
Tleaf423.ecd: leaf parameters.
Tmcrv423.ecd: microbial parameters, influencing soil decomposition, nitrogen cycle and N2O production.
Tpft423.ecd: the parameters describing the dominance of plant functional types under the control of WTD.
Tpftdecom.ecd: the parameters describing the litter decomposition rate of different plant functional types.
Trotz423.ecd: root parameters.
Tsoil423.ecd: soil density and hydrological parameters.
Tveg423_upland.ecd: vegetation parameters for peatlands.

2.3 how to run TEM

To run TEM, go to run/, make sure the intel and impi modules are loaded, and type in the command line:
mpirun -np $number_of_nodes ./xtem423e1 run.go4 tem4.log
If only running one part, then $number_of_nodes = 1. Similarly, to run 10 parts at the same time,  $number_of_nodes = 10. 
Please note that each part should be in a separate folder named as part-n/. 
To run only part-1/ but not part-0/, change to first line of para_scd.go4 to 1. 
Similarly, to run part-1 to part-10, skipping part-0, also change the first line of para_scd.go4 to 1, and $number_of_nodes = 9.
In run.go4, line 75 is how many variables you want to output, ling 76-108 is the name of these variables, and line 109-141 is the name of output files. 
If you want to output other variables, change in these three sections accordingly. The order of output variables should match the order of output file names. 
The possible output variables are in ttem423e.cpp, line 39-374. For a line like:
strcpy(predstr[I_SOLC],"SOILORGC");  
The variable name in run.go4 should be SOILORGC.

-----------------------------------------------------------------------------------------------------

3.Results: provides the key outputs of TEM simulation

3.1 Fixation: Base simulation results. 
    GFDL126fdenitmean.csv: Summed annual denitrification rate for each 0.5*0.5 grid under RCP 2.6 warming scenario. 
    GFDL126gppmean.csv: Summed annual GPP for each 0.5*0.5 grid under RCP 2.6 warming scenario. 
    GFDL126fixflmean.csv: Summed annual free-living nitrogen fixation for each 0.5*0.5 grid under RCP 2.6 warming scenario. 
    GFDL126fixsmean.csv:Summed annual symbiotic nitrogen fixation for each 0.5*0.5 grid under RCP 2.6 warming scenario. 
    GFDL126nlost.csv: Summed annual N leaching for each 0.5*0.5 grid under RCP 2.6 warming scenario.
    GFDL126soilorgc.csv:Summed annual soil carbon pool each 0.5*0.5 grid under RCP 2.6 warming scenario.
    GFDL126vegc.csv:Summed annual vegetation carbon pool each 0.5*0.5 grid under RCP 2.6 warming scenario.
    GFDL126vegnupt.csv: Summed annual vegetation N uptake for each 0.5*0.5 grid under RCP 2.6 warming scenario.
    Files staring with GFDL245 and GFDL585 represent simulation results under warming scenarios RCP 4.5  and RCP 8.5 accordingly. 

3.2 Sensitivity: Permafrost thawing simulation results.
    gfdl126n2fixflmean0.9.csv: Summed annual free-living nitrogen fixation for each 0.5*0.5 grid under RCP 2.6 warming scenario with 90% distribution. 
    gfdl126n2fixsmean0.9.csv: Summed annual symbiotic nitrogen fixation for each 0.5*0.5 grid under RCP 2.6 warming scenario with 90% distribution. 
    gfdl126n2fixflmean1.1.csv: Summed annual free-living nitrogen fixation for each 0.5*0.5 grid under RCP 2.6 warming scenario with 110% distribution. 
    gfdl126n2fixsmean1.1.csv: Summed annual symbiotic nitrogen fixation for each 0.5*0.5 grid under RCP 2.6 warming scenario with 110% distribution. 
    Files starting with gfdl245 and gfdl585represent above simulation results under warming scenarios RCP 4.5  and RCP 8.5 accordingly. 

-----------------------------------------------------------------------------------------------------

