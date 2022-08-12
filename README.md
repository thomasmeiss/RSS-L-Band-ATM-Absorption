# RSS-L-Band-ATM-Absorption

Create binary files for AO, AV, AL, TRAN, TBUP, TBDW from NCEP profiles for SMAP L2B processing

Written by Thomas Meissner + Andrew Manaster
Remote Sensing Systems
August 12, 2022

This package contains the routines for creating the atmospheric absorption ancillary fields that are ingested 
during the L2B stage of the "RSS_SMAP_Salinity-Processing-V5.0"

The routine ingests the 0.25 deg NCEP GDAS/GFS data for
1. surface temperature (2m)
2. surface relative humidity
3. surface pressure
4. surface geopotential height
5. atmospheric profile of geopotential height
6. atmospheric profile of temperature
7. atmospheric profile of relative humidity
8. atmospheric profile of cloud water mixing ratio
9. total columnar water vapor
10.total columnar cloud water

For ingestion, the NCEP fields have been extracted form the grib files (0.25 deg) and stored as raw binary with a file header containing the time stamp.

The routine calculates the columnar atmospheric absorption for oxygen water vapor and cloud liquid water absorption (AO, AV, AL) and
TRAN, TBUP, TBDW at the SMAP center frequency and Earth incidence angles 39, 40, 41 deg.
The output is stored as raw binary arrays.

Instructions:
Put all FORTRAN routines into the main directory.
Compile with MS Developer Studio.
Compile/Link Tool: compile64 (see Readme file for the "RSS-SMAP-Salinity-Processing-V5.0" package).
Put all NCEP files into the subdirectory "sample_data"
The sample data are for year=2002, month=04, day of month = 07, hour = 12Z.

Running the executable MAKE_L_BAND_ATM.exe will produce the output "atmos_ncep_2022_04_07_12z.dat" and put it into the "sample_data" directory.
This file is to be copied into the folder "ancillary_fields" of the "RSS-SMAP-Salinity-Processing-V5.0" package when running the L2B processor.
