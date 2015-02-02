import time
from datetime import datetime
import model2roms
import IOstation
import clim2bry
import DecimateGrid
import grd
import numpy as np
import atmosForcing

__author__ = 'Trond Kristiansen'
__email__ = 'trond.kristiansen@imr.no'
__created__ = datetime(2009, 1, 30)
__modified__ = datetime(2014, 12, 16)
__version__ = "1.5"
__status__ = "Development"


def myhelp():
    """
    This program is run by typing: python main.py in the command window.
    """


def showInfo(myvars, romsgridpath, climName, initName, bryName, start_year, end_year, isClimatology, useESMF, myformat):
    if isClimatology:
        print 'Conversions run for climatological months'
    else:
        print 'Conversions run from %s to year %s' % (start_year, end_year)
    print 'The following variables will be converted:'
    for myvar in myvars:
        print '---> %s' % myvar
    if (useESMF):
        print "All horisontal interpolations will be done using ESMF-ESMPy (module ESMF)"
    print "Output files are written in format: %s"%(myformat)
    print '\nOutput grid file is: %s' % romsgridpath
    print '\nInitializing'


def main():
    print '\n--------------------------\n'
    print 'Started ' + time.ctime(time.time())

    # EDIT ===================================================================
    # Set show_progress to "False" if you do not want to see the progress
    # indicator for horizontal interpolation.
    show_progress = True
    # Set compileAll to True if you want automatic re-compilation of all the
    # fortran files necessary to run soda2roms. You need to edit compile.py for this
    compileAll = False

    # Extract time-series of data for given longitude/latitude
    extractStations = False
     # Define a set of longitude/latitude positions with names to extract into
    # station files (using extractStations)
    if (extractStations):
        stationNames = ['NorthSea', 'Iceland', 'EastandWestGreenland', 'Lofoten', 'Georges Bank']
        lonlist = [2.4301, -22.6001, -47.0801, 13.3801, -67.2001]
        latlist = [54.5601, 63.7010, 60.4201, 67.5001, 41.6423]

    # Create the bry, init, and clim files for a given grid and input data
    createOceanForcing = True
    # Create atmospheric forcing for the given grid
    createAtmosForcing = False
    # Create a smaller resolution grid based on your original. Decimates every second for
    # each time run
    decimateGrid = False
    # Write ice values to file (for Arctic regions)
    writeIce = False
    # Use ESMF for the interpolation. This requires that you have ESMF and ESMPy installed (import ESMF)
    useESMF = True
    # Apply filter to smooth the 2D fields after interpolation (time consuming)
    useFilter = True

    # Format to write the ouput to: 'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT', or 'NETCDF3_CLASSIC'
    # Using NETCDF4 automatically turns on compression of files (ZLIB)
    myformat='NETCDF4'

    # Set the input data MODEL mytype
    mytype = 'SODA'
    mytype = 'SODAMONTHLY'
    mytype = 'WOAMONTHLY'
    mytype = 'NORESM'
    mytype = 'GLORYS'

    # Define what grid type you wnat to interpolate to:
    gridtype  = "NS8KM"
    #gridtype = "REGSCEN"
    #gridtype = "GREENLAND"
    #gridtype  = "KINO"

    # Define the paths to the input data
    if mytype == 'SODA':
        modelpath = "/Volumes/MacintoshHD2/Datasets/SODA/"
    if mytype == 'SODAMONTHLY':
        modelpath = "/Volumes/MacintoshHD2/Datasets/SODAMonthly/"
    if mytype == 'GLORYS':
        modelpath = "/Volumes/MacintoshHD2/Datasets/GLOBAL_REANALYSIS_PHYS_001_009/"
        modelpath = "/Users/trondkr/Projects/is4dvar/GLORYS2V3/"
    #    modelpath = "/work/users/trondk/GLORYS2V3/"
    if mytype == 'NORESM':
        modelpath = "/Users/trondkr/Projects/RegScen/NRCP45AERCN_f19_g16_CLE_01/"
        #modelpath = "/work/users/trondk/REGSCEN/NRCP45AERCN_f19_g16_CLE_01/"
        if createAtmosForcing:
            atmospath = "/Users/trondkr/Projects/RegScen/model2roms/TESTFILES/"
        
    if mytype == 'WOAMONTHLY':
        modelpath = "/Users/trondkr/Projects/is4dvar/createSSS/"

    # Define the path to the grid file
    if gridtype == "NS8KM":
        romsgridpath = "/Users/trondkr/Projects/is4dvar/Grid/nordsjoen_8km_smoothed02022015.nc"
        #romsgridpath = "/work/users/trondk/NS8km/FORCING/GRID/nordsjoen_8km_grid_hmax20m_v3.nc"

    if gridtype == "KINO":
        romsgridpath = "/work/users/trondk/KINO/GRID/kino_norseas_800m_grid.nc"
        romsgridpath="/Users/trondkr/Projects/KINO/GRID/kino_norseas_800m_grid.nc"

    if gridtype == "REGSCEN":
        romsgridpath = "/Users/trondkr/Projects/RegScen/Grid/AA_10km_grid_noest.nc"
        #romsgridpath = "/Users/trondkr/Projects/is4dvar/Grid/nordsjoen_8km_grid_hmax20m_v3.nc"
        #romsgridpath = "/work/users/trondk/REGSCEN/GRID/AA_10km_grid_noest.nc"

    if gridtype == "GREENLAND":
        romsgridpath="/Users/trondkr/Projects/RegScen/Grid/Sermilik_grid_4000m.nc"

    if mytype == 'WOAMONTHLY': isClimatology = True
    else: isClimatology = False

    # Define the period to create forcing for
    start_year  = 2009
    end_year    = 2012
    start_month = 11
    end_month   = 12

    startdate = datetime(start_year, start_month, 1)
    enddate   = datetime(end_year, end_month, 1)

    # Subset the input data. The more you subset the less memory is needed for calculations
    # and the faster the process is performed. The subset is initially performed in IOsubset.py
    if gridtype == "NS8KM":
        abbreviation = "nordsjoen_8km"
        minLat = 40
        maxLat = 70
        minLon = -30
        maxLon = 40

    if gridtype == "REGSCEN":
        abbreviation = "regscen"
        minLat = -50
        maxLat = 89.5
        minLon = -179
        maxLon = 180

    if gridtype == "GREENLAND":
        abbreviation = "greenland"
        minLat = 50
        maxLat = 89.5
        minLon = -179
        maxLon = 180

    if gridtype == "KINO":
        abbreviation = "kino"
        minLat = 30
        maxLat = 70
        minLon = -40
        maxLon = 40


    # Define what and name of variables to include in the forcing files
    # -> myvars is the name model2roms uses to identify variables
    # -> varNames is the name of the variable found in the NetCDF input files
    if mytype == 'SODA':
       myvars = ['temperature', 'salinity', 'ssh', 'uvel', 'vvel']
       varNames = ['TEMP', 'SALT', 'SSH', 'U', 'V']

    if mytype == 'SODAMONTHLY':
        myvars   = ['temperature', 'salinity', 'ssh', 'uvel', 'vvel']
        varNames = ['temp', 'salt', 'ssh', 'u', 'v']

    if mytype == 'GLORYS':
        if (writeIce):
            myvars   = ['temperature','salinity', 'ssh', 'uvel', 'vvel','uice','vice','aice','hice']
            varNames = ['votemper', 'vosaline', 'sossheig', 'vozocrtx', 'vomecrty','iicevelu', 'iicevelv', 'ileadfra', 'iicethic']
        else:
            myvars   = ['temperature', 'salinity', 'ssh', 'uvel', 'vvel']
            varNames = ['votemper', 'vosaline', 'sossheig', 'vozocrtx', 'vomecrty']

    if mytype == 'WOAMONTHLY':
        myvars   = ['temperature','salinity']
        varNames = ['t_an', 's_an']

    if mytype == 'NORESM':
        myvars   = ['temperature','salinity', 'ssh', 'uvel', 'vvel','ageice','uice','vice','aice','hice','snow_thick']
        varNames = ['templvl','salnlvl','sealv', 'uvellvl', 'vvellvl','iage', 'uvel', 'vvel', 'aice', 'hi', 'hs']

    # Define frequency of input data (5 day or 30 day average files)
    if mytype == 'SODA':
        aveDays = 5.0

    if mytype in ['SODAMONTHLY', 'GLORYS', 'NORESM','WOAMONTHLY']:
        aveDays = 30.0


    # NO EDIT BELOW =========================================================
    subset = np.zeros(4); subset[0] = minLat; subset[1] = maxLat; subset[2] = minLon; subset[3] = maxLon

    # Name of output files for CLIM, BRY, and INIT files
    climName = abbreviation + '_clim_' + str(mytype) + '_' + str(start_year) + '_to_' + str(end_year) + '.nc'
    initName = abbreviation + '_init_' + str(mytype) + '_' + str(start_year) + '_to_' + str(end_year) + '.nc'
    bryName = abbreviation + '_bry_' + str(mytype) + '_' + str(start_year) + '_to_' + str(end_year) + '.nc'
    if isClimatology is True:
        climName=abbreviation + '_' + str(mytype) + '_climatology.nc'  

    if compileAll is True:
        # Compile the Fortran 90 files to Python modules
        import compile
        compile.compileAll()

    years = [(int(startdate.year) + kk) for kk in range(1 + int(enddate.year) - int(startdate.year))]
    IDS=[]
    IDS.append(start_month)
    [IDS.append(12) for kk in range(int(enddate.year) - int(startdate.year))]
    IDS.append(end_month)

    if isClimatology==True:
        IDS=[i+1 for i in xrange(12)]
        print "Will create climatology for months: %s"%(IDS)

    # Create the grid object for the output grid
    grdROMS = grd.grdClass(romsgridpath, "ROMS", gridtype, useESMF,'ocean')
    grdROMS.vars=myvars

    if (useESMF):
        # initialize MPI
        import ESMF
        manager = ESMF.Manager(logkind = ESMF.LogKind.MULTI, debug = True)

    if createOceanForcing:

        showInfo(myvars, romsgridpath, climName, initName, bryName, start_year, end_year, isClimatology, useESMF, myformat)

        model2roms.convertMODEL2ROMS(years, IDS, climName, initName, modelpath, romsgridpath, myvars, varNames, show_progress,
                                         mytype, gridtype, subset, isClimatology, writeIce, useESMF, useFilter, myformat)

        clim2bry.writeBry(grdROMS, start_year, bryName, climName, writeIce, mytype, myformat)

    if createAtmosForcing:
        atmosForcing.createAtmosFileUV(grdROMS,modelpath,atmospath,startdate,enddate,useESMF,
            myformat,abbreviation,mytype,gridtype,show_progress)

    if decimateGrid:
        DecimateGrid.createGrid(grdROMS, '/Users/trond/Projects/arcwarm/SODA/soda2roms/imr_nordic_8km.nc', 2)

    if extractStations:
        print "Running in station mode and extracting pre-defined station locations"
        IOstation.getStationData(years, IDS, modelpath, latlist, lonlist, stationNames)

    print 'Finished ' + time.ctime(time.time())


if __name__ == "__main__":
    main()
