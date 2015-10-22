from netCDF4 import Dataset, datetime, date2num,num2date
import numpy as np
import interp2D
import interpolation as interp
import IOwrite
import os
import calendar

import grd
import barotropic
import IOinitial
import IOsubset
import datetimeFunctions

try:
    import ESMF
except ImportError:
    print "Could not find module ESMF"
    pass
__author__ = 'Trond Kristiansen'
__email__ = 'trond.kristiansen@imr.no'
__created__ = datetime(2008, 8, 15)
__modified__ = datetime(2014, 12, 1)
__version__ = "1.5"
__status__ = "Development, modified on 15.08.2008,01.10.2009,07.01.2010, 15.07.2014, 01.12.2014, 07.08.2015"

            
def VerticalInterpolation(myvar, array1, array2, grdROMS, grdMODEL):
    outINDEX_ST = (grdROMS.Nlevels, grdROMS.eta_rho, grdROMS.xi_rho)
    outINDEX_U = (grdROMS.Nlevels, grdROMS.eta_u, grdROMS.xi_u)
    outINDEX_UBAR = (grdROMS.eta_u, grdROMS.xi_u)
    outINDEX_V = (grdROMS.Nlevels, grdROMS.eta_v, grdROMS.xi_v)
    outINDEX_VBAR = (grdROMS.eta_v, grdROMS.xi_v)

    if myvar in ['salinity', 'temperature']:
        print 'Start vertical interpolation for %s (dimensions=%s x %s)' % (myvar, grdROMS.xi_rho, grdROMS.eta_rho)
        outdata = np.empty((outINDEX_ST), dtype=np.float64, order='Fortran')

        outdata = interp.interpolation.dovertinter(np.asarray(outdata, order='Fortran'),
                                                   np.asarray(array1, order='Fortran'),
                                                   np.asarray(grdROMS.h, order='Fortran'),
                                                   np.asarray(grdROMS.z_r, order='Fortran'),
                                                   np.asarray(grdMODEL.z_r, order='Fortran'),
                                                   int(grdROMS.Nlevels),
                                                   int(grdMODEL.Nlevels),
                                                   int(grdROMS.xi_rho),
                                                   int(grdROMS.eta_rho),
                                                   int(grdROMS.xi_rho),
                                                   int(grdROMS.eta_rho))

        outdata = np.ma.masked_where(abs(outdata) > 1000, outdata)

        #import plotData
        #for k in xrange(len(grdMODEL.h)-1):
        
        #    plotData.contourMap(grdROMS, grdROMS.lon_rho, grdROMS.lat_rho, np.squeeze(outdata[k,:,:]),k, myvar)


        return outdata

    if myvar == 'vvel':
        print 'Start vertical interpolation for uvel (dimensions=%s x %s)' % (grdROMS.xi_u, grdROMS.eta_u)
        outdataU = np.zeros((outINDEX_U), dtype=np.float64)
        outdataUBAR = np.zeros((outINDEX_UBAR), dtype=np.float64)

        outdataU = interp.interpolation.dovertinter(np.asarray(outdataU, order='Fortran'),
                                                    np.asarray(array1, order='Fortran'),
                                                    np.asarray(grdROMS.h, order='Fortran'),
                                                    np.asarray(grdROMS.z_r, order='Fortran'),
                                                    np.asarray(grdMODEL.z_r, order='Fortran'),
                                                    int(grdROMS.Nlevels),
                                                    int(grdMODEL.Nlevels),
                                                    int(grdROMS.xi_u),
                                                    int(grdROMS.eta_u),
                                                    int(grdROMS.xi_rho),
                                                    int(grdROMS.eta_rho))

        outdataU = np.ma.masked_where(abs(outdataU) > 1000, outdataU)

        print 'Start vertical interpolation for vvel (dimensions=%s x %s)' % (grdROMS.xi_v, grdROMS.eta_v)
        outdataV = np.zeros((outINDEX_V), dtype=np.float64)
        outdataVBAR = np.zeros((outINDEX_VBAR), dtype=np.float64)

        outdataV = interp.interpolation.dovertinter(np.asarray(outdataV, order='Fortran'),
                                                    np.asarray(array2, order='Fortran'),
                                                    np.asarray(grdROMS.h, order='Fortran'),
                                                    np.asarray(grdROMS.z_r, order='Fortran'),
                                                    np.asarray(grdMODEL.z_r, order='Fortran'),
                                                    int(grdROMS.Nlevels),
                                                    int(grdMODEL.Nlevels),
                                                    int(grdROMS.xi_v),
                                                    int(grdROMS.eta_v),
                                                    int(grdROMS.xi_rho),
                                                    int(grdROMS.eta_rho))

        outdataV = np.ma.masked_where(abs(outdataV) > 1000, outdataV)

        z_wu = np.zeros((grdROMS.Nlevels + 1, grdROMS.eta_u, grdROMS.xi_u), dtype=np.float64)
        z_wv = np.zeros((grdROMS.Nlevels + 1, grdROMS.eta_v, grdROMS.xi_v), dtype=np.float64)

        outdataUBAR = barotropic.velocity.ubar(np.asarray(outdataU, order='Fortran'),
                                               np.asarray(outdataUBAR, order='Fortran'),
                                               np.asarray(grdROMS.z_w, order='Fortran'),
                                               np.asarray(z_wu, order='Fortran'),
                                               grdROMS.Nlevels,
                                               grdROMS.xi_u,
                                               grdROMS.eta_u,
                                               grdROMS.xi_rho,
                                               grdROMS.eta_rho)
        outdataUBAR = np.ma.masked_where(abs(outdataUBAR) > 1000, outdataUBAR)

        #plotData.contourMap(grdROMS, grdROMS.lon_rho, grdROMS.lat_rho, outdataUBAR,1, "ubar")



        outdataVBAR = barotropic.velocity.vbar(np.asarray(outdataV, order='Fortran'),
                                               np.asarray(outdataVBAR, order='Fortran'),
                                               np.asarray(grdROMS.z_w, order='Fortran'),
                                               np.asarray(z_wv, order='Fortran'),
                                               grdROMS.Nlevels,
                                               grdROMS.xi_v,
                                               grdROMS.eta_v,
                                               grdROMS.xi_rho,
                                               grdROMS.eta_rho)

        #plotData.contourMap(grdROMS, grdROMS.lon_rho, grdROMS.lat_rho, outdataVBAR,1, "vbar")
        outdataVBAR = np.ma.masked_where(abs(outdataVBAR) > 1000, outdataVBAR)

        return outdataU, outdataV, outdataUBAR, outdataVBAR


def HorizontalInterpolation(myvar, grdROMS, grdMODEL, data, show_progress, useFilter):
    print 'Start %s horizontal interpolation for %s' % (grdMODEL.grdType, myvar)

    if myvar in ['temperature', 'salinity']:
        array1 = interp2D.doHorInterpolationRegularGrid(myvar, grdROMS, grdMODEL, data, show_progress, useFilter)
    if myvar in ['ssh', 'ageice', 'uice', 'vice', 'aice', 'hice', 'snow_thick']:
        array1 = interp2D.doHorInterpolationSSHRegularGrid(myvar, grdROMS, grdMODEL, data, useFilter)
    if myvar in ['uvel', 'vvel']:
        array1 = interp2D.doHorInterpolationRegularGrid(myvar, grdROMS, grdMODEL, data, show_progress, useFilter)

    return array1


def rotate(grdROMS, grdMODEL, data, u, v):
    """
    First rotate the values of U, V at rho points with the angle, and then interpolate
    the rho point values to U and V points and save the result
    """

    urot = np.zeros((int(grdMODEL.Nlevels), int(grdROMS.eta_rho), int(grdROMS.xi_rho)), np.float64)
    vrot = np.zeros((int(grdMODEL.Nlevels), int(grdROMS.eta_rho), int(grdROMS.xi_rho)), np.float64)

    urot, vrot = interp.interpolation.rotate(np.asarray(urot, order='Fortran'),
                                             np.asarray(vrot, order='Fortran'),
                                             np.asarray(u, order='Fortran'),
                                             np.asarray(v, order='Fortran'),
                                             np.asarray(grdROMS.angle, order='Fortran'),
                                             int(grdROMS.xi_rho),
                                             int(grdROMS.eta_rho),
                                             int(grdMODEL.Nlevels))
    return urot, vrot


def interpolate2UV(grdROMS, grdMODEL, urot, vrot):
    Zu = np.zeros((int(grdMODEL.Nlevels), int(grdROMS.eta_u), int(grdROMS.xi_u)), np.float64)
    Zv = np.zeros((int(grdMODEL.Nlevels), int(grdROMS.eta_v), int(grdROMS.xi_v)), np.float64)

    # Interpolate from RHO points to U and V points for velocities

    Zu = interp.interpolation.rho2u(np.asarray(Zu, order='Fortran'),
                                    np.asarray(urot, order='Fortran'),
                                    int(grdROMS.xi_rho),
                                    int(grdROMS.eta_rho),
                                    int(grdMODEL.Nlevels))


    #plotData.contourMap(grdROMS,grdMODEL,Zu[0,:,:],"1",'urot')

    Zv = interp.interpolation.rho2v(np.asarray(Zv, order='Fortran'),
                                    np.asarray(vrot, order='Fortran'),
                                    int(grdROMS.xi_rho),
                                    int(grdROMS.eta_rho),
                                    int(grdMODEL.Nlevels))


    #plotData.contourMap(grdROMS,grdMODEL,Zv[0,:,:],"1",'vrot')

    return Zu, Zv


def getTime(dataPath, indatatype, grdROMS, grdMODEL, year, month, day, mytime, firstRun):
    """
    Create a date object to keep track of Julian dates etc.
    Also create a reference date starting at 1948/01/01.
    Go here to check results:http://lena.gsfc.nasa.gov/lenaDEV/html/doy_conv.html
    """
    if indatatype == 'SODA':
        filename = getSODAfilename(year, month, day, None, dataPath)

    if indatatype == 'SODAMONTHLY':
        filename = getSODAMONTHLYfilename(year, month, day, None, dataPath)

    if indatatype == 'GLORYS':
        filename = getGLORYSfilename(year, month, day, "S", dataPath)

    if indatatype == 'WOAMONTHLY':
        filename = getWOAMONTHLYfilename(year, month, day, "temperature", dataPath)

    if indatatype == 'NORESM':
        filename = getNORESMfilename(year, month, day, "saln", dataPath)
    
    if indatatype == 'NS8KM':
        filename = getNS8KMfilename(year, month, day, "salt", dataPath)

    if indatatype == 'NS8KMZ':
        filename, readFromOneFile = getNS8KMZfilename(year, month, day, "salt", dataPath)
       
    # Now open the input file and get the time
    cdf = Dataset(filename)

    if (indatatype)=='NORESM':
        jdref = date2num(datetime(1800,1,1),cdf.variables["time"].units,calendar=cdf.variables["time"].calendar)
    elif indatatype=='NS8KMZ':
        jdref = date2num(datetime(1948,1,1),units="days since 1948-01-01 00:00:00",calendar="gregorian")
    elif (indatatype)=='GLORYS':
        jdref = date2num(datetime(1948,1,1),cdf.variables["time_counter"].units,calendar=cdf.variables["time_counter"].calendar)
    elif indatatype=='NS8KM':
        jdref = date2num(datetime(1948,1,1),cdf.variables["ocean_time"].units,calendar=cdf.variables["ocean_time"].calendar)
    else:
        jdref = date2num(datetime(1948,1,1),cdf.variables["time"].units,calendar=cdf.variables["time"].calendar)
    
    if indatatype == 'SODA':

        # Find the day and month that the SODA file respresents based on the year and ID number.
        # Each SODA file represents a 5 day average, therefore we let the date we find be the first day
        # of those 5 days. Thats the reason we subtract 4 below for day of month.
        import date

        days = 0.0; month = 1; loop = True

        while loop is True:
            d = date.NumberDaysMonth(month, year)
            if days + d < int(ID) * 5:
                days = days + d
                month += 1
            else:
                day = int(int(ID) * 5 - days)
                loop = False

        mycalendar = cdf.variables["time"].calendar
        myunits = cdf.variables["time"].units
        currentdate = datetime(year, month, day)
        jd = date2num(currentdate, myunits, calendar=mycalendar)

    if indatatype == 'SODAMONTHLY':
        # Find the day and month that the SODAMONTHLY file respresents based on the year and ID number.
        # Each SODA file represents a 1 month average.

        mycalendar = cdf.variables["time"].calendar
        myunits = cdf.variables["time"].units
        currentdate = datetime(year, month, day)
        jd = date2num(currentdate, myunits, calendar=mycalendar)


    if indatatype == 'GLORYS':
        # Find the day and month that the GLORYS file respresents based on the year and ID number.
        # Each file represents a 1 month average.
        mycalendar = cdf.variables["time_counter"].calendar
        myunits = cdf.variables["time_counter"].units
        currentdate = datetime(year, month, day)
        jd = date2num(currentdate, myunits, calendar=mycalendar)
    
    if indatatype == 'NS8KM':
        # Find the day and month that the GLORYS file respresents based on the year and ID number.
            # Each file represents a 1 month average.
        mycalendar = cdf.variables["ocean_time"].calendar
        myunits = cdf.variables["ocean_time"].units
        currentdate = datetime(year, month, day)
        jd = date2num(currentdate, myunits, calendar=mycalendar)

    if indatatype == 'NS8KMZ':
        # Find the day and month that the GLORYS file respresents based on the year and ID number.
            # Each file represents a 1 month average.
        mycalendar = "gregorian"
        refdate = datetime(1948, 1, 1)
        currentdate = datetime(year, month, day)
        myunits = cdf.variables["time"].units
        jd = date2num(currentdate, myunits, calendar="gregorian")
        print "Days:", jd, currentdate," year month day ",year, month, day

    if indatatype == 'NORESM':
        # Find the day and month that the NORESM file. We need to use the time modules from
        # netcdf4 for python as they handle calendars that are no_leap.
        # http://www.esrl.noaa.gov/psd/people/jeffrey.s.whitaker/python/netcdftime.html#datetime
        mydays = cdf.variables["time"][0]
        mycalendar = cdf.variables["time"].calendar
        myunits = cdf.variables["time"].units
        # Fake the start date first time around
       # if (firstRun):
       #     currentdate = datetime(2006,1,1)
       #     print "NOTICE!\n First datestamp in result files are hardcoded to %s"%(currentdate)
       # else:
        currentdate = num2date(mydays, units=myunits, calendar=mycalendar)
        jd = date2num(currentdate, myunits, calendar='noleap')

    grdROMS.time = (jd - jdref)
    grdROMS.reftime = jdref
    grdROMS.timeunits=myunits
    cdf.close()

    print '\nCurrent time of %s file : %s' % (indatatype, currentdate)


def getGLORYSfilename(year, month, day, myvar, dataPath):
    # Month indicates month
    # myvar:S,T,U,V
    print "MYVAR: %s - %s"%(myvar,dataPath)
    if (myvar in ['iicevelu', 'iicevelv', 'ileadfra', 'iicethic']):
        myvarPrefix = 'icemod'
        myvar = "ice"
    elif (myvar in ['sossheig']):
        if (year >= 2010 and month == 12) or ( year >= 2011 ):
            myvarPrefix = 'SSH'
        else:
            myvarPrefix = 'grid2D'
        myvar = "ssh"
    elif (myvar in ['vozocrtx', 'vomecrty']):
        myvar='u-v'
        myvarPrefix='gridUV'
    elif(myvar in ['votemper']):
        myvar='t'
        myvarPrefix='gridT'
    elif(myvar in ['vosaline']):
        myvar='s'
        myvarPrefix='gridS'
    else:
        myvarPrefix = 'grid'+str(myvar.upper())

    # GLORYS change the name in the middle of the time-series (on December 2010) and we have to
    # account for that
    if (2013 > year >= 2010 and month == 12) or ( 2013 > year >= 2011 ):
        production = "R20140520"
    elif (year >= 2013):
        production = "R20141205"
    else:
        production = "R20130808"

    if month < 10: filename = dataPath + 'dataset-global-reanalysis-phys-001-009-ran-fr-glorys2v3-monthly-' + str(
        myvar.lower()) + '/GLORYS2V3_ORCA025_' + str(year) + '0' + str(month) + '15_'+str(production)+'_' + str(myvarPrefix) + '.nc'

    if month >= 10: filename = dataPath + 'dataset-global-reanalysis-phys-001-009-ran-fr-glorys2v3-monthly-' + str(
        myvar.lower()) + '/GLORYS2V3_ORCA025_' + str(year) + str(month) + '15_'+str(production)+'_' + str(myvarPrefix) + '.nc'

    print "Filename in: %s"%(filename)

    return filename


def getNORESMfilename(year, month, day, myvar, dataPath):

    if (myvar=='grid'):
        filename = dataPath + 'GRID/NorESM.nc'
        #TODO: Fix this hardcoding of grid path
        filename = "/work/users/trondk/REGSCEN/GRID/NorESM.nc"
        filename = dataPath + 'GRID/NorESM.nc'
    else:
        if myvar in ['iage', 'uvel', 'vvel', 'aice', 'hi', 'hs']:
            filename = dataPath + 'ICE/NRCP45AERCN_f19_g16_CLE_01.cice.h.'+str(year)+'.nc'
        else:
            if (month < 10):
                filename = dataPath + 'OCN/NRCP45AERCN_f19_g16_CLE_01.micom.hm.'+str(year)+'-0'+str(ID)+'.nc'
            else:
                filename = dataPath + 'OCN/NRCP45AERCN_f19_g16_CLE_01.micom.hm.'+str(year)+'-'+str(ID)+'.nc'
 
    return filename

def getNS8KMfilename(year, month, day, myvar, dataPath):

    if (month < 10):
        mymonth='0%s'%(month)
    else:
        mymonth='%s'%(month)
    filename = dataPath + str(year)+str(mymonth)+'15_mm-IMR-MODEL-ROMS-NWS-20140430-fv02.1.nc'
      
    return filename, month

def getNS8KMZfilename(year, month, day, myvar, dataPath):

    allInOneFile = '/work/users/trondk/KINO/FORCING/1600M/northsea_8km_z_06122009_22122012.nc'
    #allInOneFile = '/Users/trondkr/Projects/KINO/northsea_8km_z_MyOceanFO_2009_2012_subset.nc'

    if os.path.exists(allInOneFile):
        readFromOneFile = True
        print "NOTE ! READING ALL MYOCEAN FORCING DATA FROM ONE FILE"
        return allInOneFile, readFromOneFile
    else:
        readFromOneFile = False
        if (month < 10):
            mymonth='0%s'%(month)
        else:
            mymonth='%s'%(month)
        filename = dataPath + str(year)+str(mymonth)+'15_mm-IMR-MODEL-ROMS-NWS-20150822-fv04.1.nc'
        return filename, readFromOneFile

def getSODAMONTHLYfilename(year, month, day, myvar, dataPath):

    if month < 10: filename = dataPath + 'SODA_2.0.2_' + str(year) + '0' + str(month) + '.cdf'
    if month >= 10: filename = dataPath + 'SODA_2.0.2_' + str(year) + str(month) + '.cdf'

    return filename

def getSODAfilename(year, month, day, myvar, dataPath):

    file = "SODA_2.0.2_" + str(year) + "_" + str(month) + ".cdf"

    return dataPath + file

def getWOAMONTHLYfilename(year, month, day, myvar, dataPath):

    if myvar == "temperature":
        filename = dataPath + 'temperature_monthly_1deg.nc'
    elif myvar == "salinity":
        filename = dataPath + 'salinity_monthly_1deg.nc'
    else:
        print "Could not find any input files in folder: %s"%(datapath)

    return filename


def get3Ddata(grdROMS, grdMODEL, myvar, indatatype, year, month, day, varNames, dataPath):

    if myvar == 'temperature': varN = 0;
    if myvar == 'salinity': varN = 1;
    if myvar == 'uvel': varN = 3;
    if myvar == 'vvel': varN = 4;

    # The variable splitExtract is defined in IOsubset.py and depends on the orientation
    # and indatatype of grid (-180-180 or 0-360). Assumes regular grid.
    if  grdMODEL.useESMF:
        if indatatype == "SODA":
            filename = getSODAfilename(year, month, day, None, dataPath)
            cdf = Dataset(filename)
            data = cdf.variables[varNames[varN]][0,:,:,:]

        if indatatype == "SODAMONTHLY":
            filename = getSODAMONTHLYfilename(year, month, day, None, dataPath)
            cdf = Dataset(filename)
            data = cdf.variables[str(varNames[varN])][:,:,:]

        if indatatype == "WOAMONTHLY":
            filename = getWOAMONTHLYfilename(year, month, day, myvar, dataPath)
            cdf = Dataset(filename)
            data = cdf.variables[str(varNames[varN])][month - 1, :,:,:]

        if indatatype == "NORESM":
            cdf = Dataset(getNORESMfilename(year, month, day, varNames[varN], dataPath))
            myunits = cdf.variables[str(varNames[varN])].units
            data = np.squeeze(cdf.variables[str(varNames[varN])][0,:,:,:])
            data=np.where(data.mask,grdROMS.fill_value,data)
           # data = np.where(abs(data)>=32768 , grdROMS.fill_value, data)

            print "Data range", np.min(data),np.max(data)

        if indatatype == "NS8KMZ":
            filename, readFromOneFile = getNS8KMZfilename(year, month, day, varNames[varN], dataPath)
            cdf = Dataset(filename)
            print "Reading from one file %s"%(readFromOneFile)

            myunits = cdf.variables[str(varNames[varN])].units
            if (readFromOneFile):
                jdref = date2num(datetime(1948,1,1),cdf.variables["time"].units,calendar=cdf.variables["time"].calendar)
                currentdate = datetime(year, month, day, 12)
                jd = date2num(currentdate, cdf.variables["time"].units, calendar="gregorian")
                print "Days:", jd, currentdate
  
                timesteps = (cdf.variables["time"][:]).tolist()
                timeindex = timesteps.index(jd)
                data = np.squeeze(cdf.variables[str(varNames[varN])][timeindex,:,:,:])
            else:
                data = np.squeeze(cdf.variables[str(varNames[varN])][0,:,:,:])
                print "Range of data",varN, np.min(data),np.max(data)
            data=np.where(data.mask,grdROMS.fill_value,data)
     
        if indatatype == "GLORYS":
            cdf = Dataset(getGLORYSfilename(year, month, day, varNames[varN], dataPath))
            myunits = cdf.variables[str(varNames[varN])].units
            data = np.squeeze(cdf.variables[str(varNames[varN])][0,:,:,:])
            data=np.where(data.mask,grdROMS.fill_value,data)

        cdf.close()
    else:
        if grdMODEL.splitExtract is True:
            if indatatype == "SODA":
                filename = getSODAfilename(year, month, day,  None, dataPath)
                cdf = Dataset(filename)

                data1 = cdf.variables[varNames[varN]][0, :,
                        int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                        int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]
                data2 = cdf.variables[varNames[varN]][0, :,
                        int(grdMODEL.indices[1, 2]):int(grdMODEL.indices[1, 3]),
                        int(grdMODEL.indices[1, 0]):int(grdMODEL.indices[1, 1])]

            if indatatype == "SODAMONTHLY":
                filename = getSODAMONTHLYfilename(year, month, day, None, dataPath)
                cdf = Dataset(filename)

                data1 = cdf.variables[str(varNames[varN])][:,
                        int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                        int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]
                data2 = cdf.variables[str(varNames[varN])][:,
                        int(grdMODEL.indices[1, 2]):int(grdMODEL.indices[1, 3]),
                        int(grdMODEL.indices[1, 0]):int(grdMODEL.indices[1, 1])]

            if indatatype in ["WOAMONTHLY"]:

                filename = getWOAMONTHLYfilename(year, month, day, myvar, dataPath)
                cdf = Dataset(filename)

                data1 = cdf.variables[str(varNames[varN])][month - 1, :,
                        int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                        int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]
                data2 = cdf.variables[str(varNames[varN])][month - 1, :,
                        int(grdMODEL.indices[1, 2]):int(grdMODEL.indices[1, 3]),
                        int(grdMODEL.indices[1, 0]):int(grdMODEL.indices[1, 1])]

            if indatatype == "GLORYS":
                cdf = Dataset(getGLORYSfilename(year, month, day, varNames[varN], dataPath))
                myunits = cdf.variables[str(varNames[varN])].units

                data1 = np.squeeze(cdf.variables[str(varNames[varN])][0, :,
                                   int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                                   int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])])
                data2 = np.squeeze(cdf.variables[str(varNames[varN])][0, :,
                                   int(grdMODEL.indices[1, 2]):int(grdMODEL.indices[1, 3]),
                                   int(grdMODEL.indices[1, 0]):int(grdMODEL.indices[1, 1])])

            if indatatype == "NORESM":
                cdf = Dataset(getNORESMfilename(year, month, day, varNames[varN], dataPath))
                myunits = cdf.variables[str(varNames[varN])].units
                data1 = np.squeeze(cdf.variables[str(varNames[varN])][0, :,
                                   int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                                   int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])])
                data2 = np.squeeze(cdf.variables[str(varNames[varN])][0, :,
                                   int(grdMODEL.indices[1, 2]):int(grdMODEL.indices[1, 3]),
                                   int(grdMODEL.indices[1, 0]):int(grdMODEL.indices[1, 1])])

            cdf.close()

            data = np.concatenate((data1, data2), axis=2)

        else:
            if indatatype == "SODA":
                filename = getSODAfilename(year, month, day, None, dataPath)
                cdf = Dataset(filename)

                data = cdf.variables[str(varNames[varN])][0, :,
                       int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                       int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]

            if indatatype == "SODAMONTHLY":
                filename = getSODAMONTHLYfilename(year, month, day, None, dataPath)
                cdf = Dataset(filename)

                data = cdf.variables[str(varNames[varN])][:,
                       int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                       int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]

            if indatatype == "WOAMONTHLY":
                filename = getWOAMONTHLYfilename(year, month, day, myvar, dataPath)
                cdf = Dataset(filename)

                data = cdf.variables[str(varNames[varN])][month, :,
                       int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                       int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]

            if indatatype == "GLORYS":
                if myvar == 'temperature':
                    cdf = Dataset(getGLORYSfilename(year, month, day, 'T', dataPath))

                    myunits = cdf.variables[str(varNames[varN])].units
                if myvar == 'salinity': cdf = Dataset(getGLORYSfilename(year, month, day, 'S', dataPath))
                if myvar == 'uvel': cdf = Dataset(getGLORYSfilename(year, month, day, 'U', dataPath))
                if myvar == 'vvel': cdf = Dataset(getGLORYSfilename(year, month, day, 'V', dataPath))

                data = np.squeeze(cdf.variables[str(varNames[varN])][0, :,
                                  int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                                  int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])])

            if indatatype == "NORESM":
                cdf = Dataset(getNORESMfilename(year, month, day, varNames[varN], dataPath))
                myunits = cdf.variables[str(varNames[varN])].units

                data = np.squeeze(cdf.variables[str(varNames[varN])][0, :,
                                  int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                                  int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])])
            cdf.close()

    if myvar == 'temperature' and indatatype in ["NS8KMZ", "GLORYS", "NORESM"]:

        if myunits == "degree_Kelvin" or myunits == "K":
            if indatatype in ["GLORYS"]:
                data = np.where(data <= -32.767, grdROMS.fill_value, data)
            data = data - 273.15

            #  if time == 0 and myvar == myvars[0]:
            #      tmp = np.squeeze(data[0, :, :])
            #grdMODEL.mask = np.zeros(grdMODEL.lon.shape, dtype=np.float64)
            #grdMODEL.mask[:, :] = np.where(tmp == grdROMS.fill_value, 1, 0)

    if indatatype == "GLORYS":
        data = np.where(data <= -32.767, grdROMS.fill_value, data)
        data = np.ma.masked_where(data <= grdROMS.fill_value, data)


    if __debug__:
        print "Data range of %s just after extracting from netcdf file: %s - %s" % (str(varNames[varN]),
                                                                                    data.min(), data.max())

    return data

def get2Ddata(grdROMS, grdMODEL, myvar, indatatype, year, month, day, varNames, dataPath):

    indexROMS_SSH = (grdROMS.eta_rho, grdROMS.xi_rho)

    if indatatype == "NORESM":
        if myvar == 'ageice': varN = 5;  SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'uice':   varN = 6;  SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'vice':   varN = 7;  SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'aice':   varN = 8;  SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'hice':   varN = 9;  SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'snow_thick': varN = 10;  SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)

    if indatatype == "GLORYS":
        if myvar == 'uice':   varN = 5;  SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'vice':   varN = 6;  SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'aice':   varN = 7;  SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)
        if myvar == 'hice':   varN = 8;  SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)

    if myvar == 'ssh': varN = 2;  SSHdata = np.zeros((indexROMS_SSH), dtype=np.float64)

    if  grdMODEL.useESMF:
        if indatatype == "SODA":
            filename = getSODAfilename(year, month, day, None, dataPath)
            cdf = Dataset(filename)
            data = cdf.variables[varNames[varN]][0,:,:]

        if indatatype == "SODAMONTHLY":
            filename = getSODAMONTHLYfilename(year, month, day, None, dataPath)
            cdf = Dataset(filename)
            data = cdf.variables[str(varNames[varN])][:,:]

        if indatatype == "WOAMONTHLY":
            filename = getWOAMONTHLYfilename(year, month, day, myvar, dataPath)
            cdf = Dataset(filename)
            data = cdf.variables[str(varNames[varN])][month - 1, :,:]

        if indatatype == "NORESM":
            cdf = Dataset(getNORESMfilename(year, month, day, varNames[varN], dataPath))
            #myunits = cdf.variables[str(varNames[varN])].units
            # For NORESM data are 12 months of data stored in ice files. Use ID as month indicator to get data.
            if myvar in ['ageice','uice','vice','aice','hice','snow_thick']:
                data = np.squeeze(cdf.variables[str(varNames[varN])][month-1, :,:])
            else:
                data = np.squeeze(cdf.variables[str(varNames[varN])][0, :,:])
            data=np.where(data.mask,grdROMS.fill_value,data)

        if indatatype == "GLORYS":
            cdf = Dataset(getGLORYSfilename(year, month, day, varNames[varN], dataPath))
            data = np.squeeze(cdf.variables[str(varNames[varN])][0, :,:])
            data=np.where(data.mask,grdROMS.fill_value,data)

        if indatatype == "NS8KMZ":
            filename, readFromOneFile = getNS8KMZfilename(year, month, day, varNames[varN], dataPath)
            cdf = Dataset(filename)
            print "Reading from file", filename
            if (readFromOneFile):
                jdref = date2num(datetime(1948,1,1),cdf.variables["time"].units,calendar=cdf.variables["time"].calendar)
                currentdate = datetime(year, month, day, 12)
                jd = date2num(currentdate, cdf.variables["time"].units, calendar="gregorian")
              
                timesteps = (cdf.variables["time"][:]).tolist()
                timeindex = timesteps.index(jd)
                data = np.squeeze(cdf.variables[str(varNames[varN])][timeindex,:,:])
            else:
                data = np.squeeze(cdf.variables[str(varNames[varN])][0,:,:])
                
                
            data=np.where(data.mask,grdROMS.fill_value,data)
            
            print "Extracted raw data: %s min: %s max: %s"%(myvar,np.min(data),np.max(data))
        cdf.close()
    else:
        if grdMODEL.splitExtract is True:
            if indatatype == "SODA":
                filename = getSODAfilename(year, month, day, None, dataPath)
                cdf = Dataset(filename)

                data1 = cdf.variables[varNames[varN]][0,
                        int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                        int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]
                data2 = cdf.variables[varNames[varN]][0,
                        int(grdMODEL.indices[1, 2]):int(grdMODEL.indices[1, 3]),
                        int(grdMODEL.indices[1, 0]):int(grdMODEL.indices[1, 1])]

            if indatatype == "SODAMONTHLY":
                filename = getSODAMONTHLYfilename(year, month, day, None, dataPath)
                cdf = Dataset(filename)

                data1 = cdf.variables[str(varNames[varN])][
                        int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                        int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]
                data2 = cdf.variables[str(varNames[varN])][
                        int(grdMODEL.indices[1, 2]):int(grdMODEL.indices[1, 3]),
                        int(grdMODEL.indices[1, 0]):int(grdMODEL.indices[1, 1])]

            if indatatype == "GLORYS":
                cdf = Dataset(getGLORYSfilename(year, month, day, '2D', dataPath))

                data1 = np.squeeze(cdf.variables[str(varNames[varN])][0,
                                   int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                                   int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])])
                data2 = np.squeeze(cdf.variables[str(varNames[varN])][0,
                                   int(grdMODEL.indices[1, 2]):int(grdMODEL.indices[1, 3]),
                                   int(grdMODEL.indices[1, 0]):int(grdMODEL.indices[1, 1])])

            if indatatype == "NORESM":
                cdf = Dataset(getNORESMfilename(year, month, day, varNames[varN], dataPath))

                if myvar in ['ageice','uice','vice','aice','hice','snow_thick']:
                    data1 = np.squeeze(cdf.variables[str(varNames[varN])][month-1,
                                   int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                                   int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])])
                    data2 = np.squeeze(cdf.variables[str(varNames[varN])][month-1,
                                   int(grdMODEL.indices[1, 2]):int(grdMODEL.indices[1, 3]),
                                   int(grdMODEL.indices[1, 0]):int(grdMODEL.indices[1, 1])])
                else:
                    data1 = np.squeeze(cdf.variables[str(varNames[varN])][0,
                                   int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                                   int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])])
                    data2 = np.squeeze(cdf.variables[str(varNames[varN])][0,
                                   int(grdMODEL.indices[1, 2]):int(grdMODEL.indices[1, 3]),
                                   int(grdMODEL.indices[1, 0]):int(grdMODEL.indices[1, 1])])

            cdf.close()
            data = np.concatenate((data1, data2), axis=1)

        else:
            if indatatype == "SODA":
                filename = getSODAfilename(year, month, day, None, dataPath)
                cdf = Dataset(filename)

                data = cdf.variables[str(varNames[varN])][0,
                       int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                       int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]

            if indatatype == "SODAMONTHLY":
                filename = getSODAMONTHLYfilename(year, month, day, None, dataPath)
                cdf = Dataset(filename)

                data = cdf.variables[str(varNames[varN])][
                       int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                       int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])]

            if indatatype == "GLORYS":
                cdf = Dataset(getGLORYSfilename(year, month, day, varNames[varN], dataPath))

                data = np.squeeze(cdf.variables[str(varNames[varN])][0,
                                  int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                                  int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])])

            if indatatype == "NORESM":
                cdf = Dataset(getNORESMfilename(year, month, day, varNames[varN], dataPath))

                if myvar in ['ageice','uice','vice','aice','hice','snow_thick']:
                    data = np.squeeze(cdf.variables[str(varNames[varN])][month-1,
                                  int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                                  int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])])
                else:
                    data = np.squeeze(cdf.variables[str(varNames[varN])][0,
                                  int(grdMODEL.indices[0, 2]):int(grdMODEL.indices[0, 3]),
                                  int(grdMODEL.indices[0, 0]):int(grdMODEL.indices[0, 1])])

            cdf.close()
            # data = np.where(abs(data) == grdMODEL.fill_value, grdROMS.fill_value, data)
            #data = np.ma.where(abs(data) > 1000, grdROMS.fill_value, data)

        if indatatype == "GLORYS":
            data = np.where(data <= -32.767, grdROMS.fill_value, data)
            data = np.ma.masked_where(data <= grdROMS.fill_value, data)

        if __debug__:
            print "Data range of %s just after extracting from netcdf file: %s - %s" % (str(varNames[varN]),
                                                                                        data.min(), data.max())

    return data


def convertMODEL2ROMS(years, startdate, enddate, timeFrequencyOfInputData, climName, initName, dataPath, romsgridpath, myvars, varNames, show_progress, indatatype, gridtype,
                      isClimatology, writeIce, useESMF, useFilter, myformat, subsetIndata, outgrid, subset):
    if useESMF: 
        print "Turning on debugging for ESMF"
        ESMF.Manager(logkind=ESMF.LogKind.MULTI, debug=True)
    
    # First opening of input file is just for initialization of grid
    if indatatype == 'SODA':
        fileNameIn = getSODAfilename(startdate.year, startdate.month, startdate.day, "salinity", dataPath)
    if indatatype == 'SODAMONTHLY':
        fileNameIn = getSODAfilename(startdate.year, startdate.month, startdate.day, "salinity", dataPath)
    if indatatype == 'NORESM':
        fileNameIn = getNORESMfilename(startdate.year, startdate.month, startdate.day, "grid", dataPath)
    if indatatype == 'WOAMONTHLY':
        fileNameIn = getWOAMONTHLYfilename(startdate.year, startdate.month, startdate.day, "temperature", dataPath)
    if indatatype == 'GLORYS':
        fileNameIn = getGLORYSfilename(startdate.year, startdate.month, startdate.day, "S", dataPath)
    if indatatype == 'GLORYS':
        fileNameIn = getGLORYSfilename(startdate.year, startdate.month, startdate.day, "S", dataPath)
    if indatatype == 'NS8KM':
        fileNameIn = getNS8KMfilename(startdate.year, startdate.month, startdate.day, "S", dataPath)
    if indatatype == 'NS8KMZ':
        fileNameIn, readFromOneFile = getNS8KMZfilename(startdate.year, startdate.month, startdate.day, "S", dataPath)
     
    # First time in loop, get the essential old grid information
    # MODEL data already at Z-levels. No need to interpolate to fixed depths,
    # but we use the one we have

    grdMODEL = grd.grdClass(fileNameIn, indatatype, indatatype, useESMF,"ocean", outgrid)
    grdROMS = grd.grdClass(romsgridpath, "ROMS", gridtype, useESMF,"ocean", outgrid)
    grdROMS.myvars = myvars
    if (useESMF):
        print "\nCreating the interpolation weights and indexes using ESMF (this may take some time....):"

        print "  -> regridSrc2Dst at RHO points"
        grdMODEL.fieldSrc = ESMF.Field(grdMODEL.esmfgrid, "fieldSrc", staggerloc=ESMF.StaggerLoc.CENTER)
        grdMODEL.fieldDst_rho = ESMF.Field(grdROMS.esmfgrid, "fieldDst", staggerloc=ESMF.StaggerLoc.CENTER)
        grdMODEL.regridSrc2Dst_rho = ESMF.Regrid(grdMODEL.fieldSrc, grdMODEL.fieldDst_rho, regrid_method=ESMF.RegridMethod.BILINEAR, unmapped_action=ESMF.UnmappedAction.IGNORE)
      
        print "  -> regridSrc2Dst at U points"
        grdMODEL.fieldSrc = ESMF.Field(grdMODEL.esmfgrid, "fieldSrc", staggerloc=ESMF.StaggerLoc.CENTER)
        grdMODEL.fieldDst_u = ESMF.Field(grdROMS.esmfgrid_u, "fieldDst_u", staggerloc=ESMF.StaggerLoc.CENTER)
        grdMODEL.regridSrc2Dst_u = ESMF.Regrid(grdMODEL.fieldSrc, grdMODEL.fieldDst_u, regrid_method=ESMF.RegridMethod.BILINEAR, unmapped_action=ESMF.UnmappedAction.IGNORE)

        print "  -> regridSrc2Dst at V points"
        grdMODEL.fieldSrc = ESMF.Field(grdMODEL.esmfgrid, "fieldSrc", staggerloc=ESMF.StaggerLoc.CENTER)
        grdMODEL.fieldDst_v = ESMF.Field(grdROMS.esmfgrid_v, "fieldDst_v", staggerloc=ESMF.StaggerLoc.CENTER)
        grdMODEL.regridSrc2Dst_v = ESMF.Regrid(grdMODEL.fieldSrc, grdMODEL.fieldDst_v, regrid_method=ESMF.RegridMethod.BILINEAR, unmapped_action=ESMF.UnmappedAction.IGNORE)

    # Now we want to subset the data to avoid storing more information than we need.
    # We do this by finding the indices of maximum and minimum latitude and longitude in the matrixes
    if subsetIndata:
        IOsubset.findSubsetIndices(grdMODEL, min_lat=subset[0], max_lat=subset[1], min_lon=subset[2], max_lon=subset[3])

    print 'Initializing done'
    print '\n--------------------------'
    time = 0; firstRun = True

    for year in years:
        months = datetimeFunctions.createListOfMonths(year,startdate,enddate,isClimatology)
        
        for month in months:
            days = datetimeFunctions.createListOfDays(year,month,startdate,enddate,isClimatology,timeFrequencyOfInputData)

            for day in days:
                # Get the current date for given timestep 
                getTime(dataPath, indatatype, grdROMS, grdMODEL, year, month, day, time, firstRun)
               
                # Each MODEL file consist only of one time step. Get the subset data selected, and
                # store that time step in a new array:

                if firstRun is True:
                    print "NOTE! Make sure that these two arrays are in sequential order:"
                    print "myvars:     %s" % (myvars)
                    print "varnames: %s\n" % (varNames)
                    firstRun = False

                    if subsetIndata:
                        # The first iteration we want to organize the subset indices we want to extract
                        # from the input data to get the interpolation correct and to function fast
                        IOsubset.organizeSplit(grdMODEL, grdROMS)

                for myvar in myvars:

                    if myvar in ['temperature', 'salinity', 'uvel', 'vvel']:
                        data = get3Ddata(grdROMS, grdMODEL, myvar, indatatype, year, month, day, varNames, dataPath)

                    if myvar in ['ssh', 'ageice', 'uice', 'vice', 'aice', 'hice','snow_thick']:
                        data = get2Ddata(grdROMS, grdMODEL, myvar, indatatype, year, month, day, varNames, dataPath)

                    # Take the input data and horizontally interpolate to your grid
                    array1 = HorizontalInterpolation(myvar, grdROMS, grdMODEL, data, show_progress, useFilter)
                    if myvar in ['temperature', 'salinity']:
                        STdata = VerticalInterpolation(myvar, array1, array1, grdROMS, grdMODEL)
                        print "Data range of %s after interpolation: %3.3f to %3.3f" % (
                            myvar, STdata.min(), STdata.max())

                        for dd in xrange(len(STdata[:,0,0])):
                            STdata[dd,:,:] = np.where(grdROMS.mask_rho == 0, grdROMS.fill_value, STdata[dd,:,:])

                        STdata = np.where(abs(STdata) > 1000, grdROMS.fill_value, STdata)

                        IOwrite.writeClimFile(grdROMS, time, climName, myvar, isClimatology, writeIce, indatatype, myformat, STdata)
                        if time == grdROMS.initTime and grdROMS.write_init is True:
                            IOinitial.createInitFile(grdROMS, time, initName, myvar, writeIce, indatatype, myformat, STdata)

                    if myvar in ['ssh', 'ageice', 'aice', 'hice', 'snow_thick']:
                        SSHdata = array1[0, :, :]

                        SSHdata = np.where(grdROMS.mask_rho == 0, grdROMS.fill_value, SSHdata)
                        SSHdata = np.where(abs(SSHdata) > 100, grdROMS.fill_value, SSHdata)
                        SSHdata = np.where(abs(SSHdata) == 0, grdROMS.fill_value, SSHdata)

                        # Specific for ROMs. We set 0 where we should have fillvalue for ice otherwise ROMS blows up.
                        SSHdata = np.where(abs(SSHdata) == grdROMS.fill_value, 0, SSHdata)
                       # SSHdata = np.ma.masked_where(abs(SSHdata) > 100, SSHdata)

                       # print "Data range of %s after interpolation: %3.3f to %3.3f" % (
                       #     myvar, SSHdata.min(), SSHdata.max())

                        IOwrite.writeClimFile(grdROMS, time, climName, myvar, isClimatology, writeIce, indatatype,myformat, SSHdata)
                        if time == grdROMS.initTime:
                            IOinitial.createInitFile(grdROMS, time, initName, myvar, writeIce, indatatype, myformat,  SSHdata)

                    # The following are special routines used to calculate the u and v velocity
                    # of ice based on the transport, which is divided by snow and ice thickenss
                    # and then multiplied by grid size in dx or dy direction (opposite of transport).
                    if myvar in ['uice', 'vice']:
                        SSHdata = array1[0, :, :]

                        if  myvar=="uice":mymask=grdROMS.mask_u
                        if  myvar=="vice":mymask=grdROMS.mask_v

                        SSHdata = np.where(mymask == 0, grdROMS.fill_value, SSHdata)
                        SSHdata = np.where(abs(SSHdata) > 100, grdROMS.fill_value, SSHdata)
                        SSHdata = np.where(abs(SSHdata) == 0, grdROMS.fill_value, SSHdata)
                        SSHdata = np.where(abs(SSHdata) == grdROMS.fill_value, 0, SSHdata)

                        #SSHdata = np.ma.masked_where(abs(SSHdata) > 1000, SSHdata)

                        print "Data range of %s after interpolation: %3.3f to %3.3f" % (myvar, SSHdata.min(), SSHdata.max())


                        IOwrite.writeClimFile(grdROMS, time, climName, myvar, isClimatology, writeIce, indatatype, myformat, SSHdata)

                        if time == grdROMS.initTime:
                            if myvar == 'uice':
                                IOinitial.createInitFile(grdROMS, time, initName, 'uice', writeIce, indatatype, myformat,  SSHdata)
                            if myvar == 'vice':
                                IOinitial.createInitFile(grdROMS, time, initName, 'vice', writeIce, indatatype, myformat, SSHdata)

                    if myvar == 'uvel':
                        array2 = array1

                    if myvar == 'vvel':
                        indexROMS_UBAR = (grdROMS.eta_u, grdROMS.xi_u)
                        indexROMS_VBAR = (grdROMS.eta_v, grdROMS.xi_v)
                        UBARdata = np.zeros((indexROMS_UBAR), dtype=np.float64)
                        VBARdata = np.zeros((indexROMS_VBAR), dtype=np.float64)

                        urot, vrot = rotate(grdROMS, grdMODEL, data, array2, array1)

                        u, v = interpolate2UV(grdROMS, grdMODEL, urot, vrot)

                        Udata, Vdata, UBARdata, VBARdata = VerticalInterpolation(myvar, u, v, grdROMS, grdMODEL)

                    if myvar == 'vvel':
                     #   print "Data range of U after interpolation: %3.3f to %3.3f - V after scaling: %3.3f to %3.3f" % (
                      #      Udata.min(), Udata.max(), Vdata.min(), Vdata.max())

                        Udata = np.where(grdROMS.mask_u == 0, grdROMS.fill_value, Udata)
                        Udata = np.where(abs(Udata) > 1000, grdROMS.fill_value, Udata)
                        Vdata = np.where(grdROMS.mask_v == 0, grdROMS.fill_value, Vdata)
                        Vdata = np.where(abs(Vdata) > 1000, grdROMS.fill_value, Vdata)
                        UBARdata = np.where(grdROMS.mask_u == 0, grdROMS.fill_value, UBARdata)
                        UBARdata = np.where(abs(UBARdata) > 1000, grdROMS.fill_value, UBARdata)
                        VBARdata = np.where(grdROMS.mask_v == 0, grdROMS.fill_value, VBARdata)
                        VBARdata = np.where(abs(VBARdata) > 1000, grdROMS.fill_value, VBARdata)

                        IOwrite.writeClimFile(grdROMS, time, climName, myvar, isClimatology, writeIce, indatatype, myformat, Udata, Vdata,
                                              UBARdata, VBARdata)
                        if time == grdROMS.initTime:
                            # We print time=initTime to init file so that we have values for ubar and vbar (not present at time=1)
                            IOinitial.createInitFile(grdROMS, time, initName, myvar, writeIce, indatatype, myformat, Udata, Vdata, UBARdata, VBARdata)

                time += 1
