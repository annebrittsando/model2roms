
from netCDF4 import Dataset
from netCDF4 import num2date
import numpy as np
import time
import os
import shutil

def createGrid(grdROMS,infilename,outfilename,decimate):
        
    shutil.copy2(infilename, outfilename)
    
    f1 = Dataset(outfilename, mode='w', format='NETCDF4')
    f1.description="This is a grid file for ROMS - KINO project"
    f1.history = 'Created (decimated) from IMR 800M grid' + time.ctime(time.time())
    f1.source = 'Trond Kristiansen (trond.kristiansen@imr.no)'
    f1.type='NetCDF4 classic created using MODEL2ROMS - https://github.com/trondkr/model2roms'
    
    if not int(grdROMS.xi_rho) % 2 == 0: 
        deltaXI=1
    else:
        deltaXI=0
    if not int(grdROMS.eta_rho) % 2 == 0: 
        deltaETA=1
        startindex=0;endindex=-2
    else:
        deltaETA=0
        startindex=0;endindex=-1

    print 'Old dimensions were    : %ix%i'%(int(grdROMS.xi_rho)-deltaXI,int(grdROMS.eta_rho)-deltaETA)
    print 'Old dimensions were    : %ix%i'%(int(grdROMS.xi_u)-deltaXI,int(grdROMS.eta_u)-deltaETA)
    print 'Old dimensions were    : %ix%i'%(int(grdROMS.xi_v)-deltaXI,int(grdROMS.eta_v)-deltaETA)
    print 'Old dimensions were    : %ix%i'%(int(grdROMS.xi_psi)-deltaXI,int(grdROMS.eta_psi)-deltaETA)
     
    print 'New dimensions will be : %ix%i'%(((int(grdROMS.xi_rho)-deltaXI)/decimate),((int(grdROMS.eta_rho)-deltaETA)/decimate))
    print 'New dimensions will be : %ix%i'%(((int(grdROMS.xi_u)-deltaXI)/decimate),((int(grdROMS.eta_u)-deltaETA)/decimate))
    print 'New dimensions will be : %ix%i'%(((int(grdROMS.xi_v)-deltaXI)/decimate),((int(grdROMS.eta_v)-deltaETA)/decimate))
    print 'New dimensions will be : %ix%i'%(((int(grdROMS.xi_psi)-deltaXI)/decimate),((int(grdROMS.eta_psi)-deltaETA)/decimate))
    
    # Define dimensions
    xi_rho=(int(grdROMS.xi_rho)-deltaXI)/decimate
    xi_vert=((int(grdROMS.xi_rho)-deltaXI)/decimate) + 1
    eta_rho=(int(grdROMS.eta_rho)-deltaETA)/decimate
    eta_vert=((int(grdROMS.eta_rho)-deltaETA)/decimate) + 1 
    xi_u=(int(grdROMS.xi_u)-deltaXI)/decimate
    eta_u=(int(grdROMS.eta_u)-deltaETA)/decimate
    xi_v=(int(grdROMS.xi_v)-deltaXI)/decimate
    eta_v=(int(grdROMS.eta_v)-deltaETA)/decimate

    xi_psi=(int(grdROMS.xi_psi)-deltaXI)/decimate
    eta_psi=(int(grdROMS.eta_psi)-deltaETA)/decimate
    s_rho=int(len(grdROMS.s_rho))
    s_w=int(len(grdROMS.s_w))

    f1.createDimension('xi_rho',  xi_rho)
    f1.createDimension('eta_rho', eta_rho)
    f1.createDimension('xi_u',    xi_u)
    f1.createDimension('eta_u',   eta_u)
    f1.createDimension('xi_v',    xi_v)
    f1.createDimension('eta_v',   eta_v)
    f1.createDimension('xi_psi',  xi_psi)
    f1.createDimension('eta_psi', eta_psi)
    f1.createDimension('xi_vert', xi_vert)
    f1.createDimension('eta_vert', eta_vert)
    f1.createDimension('s_rho',   s_rho)
    f1.createDimension('s_w',     s_w)

    vnc = f1.createVariable('lon_rho', 'd', ('eta_rho','xi_rho',),zlib=True)
    vnc.long_name = 'Longitude at RHO points'
    vnc.units = 'degrees east'
    vnc[:,:] = grdROMS.lon_rho[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc = f1.createVariable('lat_rho', 'd', ('eta_rho','xi_rho',),zlib=True)
    vnc.long_name = 'Latitude at RHO points'
    vnc.units = 'degrees north'
    vnc[:,:] = grdROMS.lat_rho[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc = f1.createVariable('lon_u', 'd', ('eta_u','xi_u',),zlib=True)
    vnc.long_name = 'Longitude at U points'
    vnc.units = 'degrees east'
    vnc[:,:] = grdROMS.lon_u[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc = f1.createVariable('lat_u', 'd', ('eta_u','xi_u',),zlib=True)
    vnc.long_name = 'Latitude at U points'
    vnc.units = 'degrees north'
    vnc[:,:] = grdROMS.lat_u[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc = f1.createVariable('lon_v', 'd', ('eta_v','xi_v',),zlib=True)
    vnc.long_name = 'Longitude at V points'
    vnc.units = 'degrees east'
    vnc[:,:] = grdROMS.lon_v[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc = f1.createVariable('lat_v', 'd', ('eta_v','xi_v',),zlib=True)
    vnc.long_name = 'Latitude at V points'
    vnc.units = 'degrees north'
    vnc[:,:] = grdROMS.lat_v[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc = f1.createVariable('lat_psi', 'd', ('eta_psi','xi_psi',),zlib=True)
    vnc.long_name = 'Latitude at PSI points'
    vnc.units = 'degrees north'
    vnc[:,:] = grdROMS.lat_psi[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc = f1.createVariable('lon_psi', 'd', ('eta_psi','xi_psi',),zlib=True)
    vnc.long_name = 'Longitude at PSI points'
    vnc.units = 'degrees east'
    vnc[:,:] = grdROMS.lon_psi[startindex:endindex:decimate,startindex:endindex:decimate]
   
    vnc = f1.createVariable('lon_vert', 'd', ('eta_vert','xi_vert',),zlib=True)
    vnc.long_name = 'Longitude at vertices points'
    vnc.units = 'degrees east'
    vnc[:,:] = grdROMS.lon_vert[::decimate,::decimate]
    
    vnc = f1.createVariable('lat_vert', 'd', ('eta_vert','xi_vert',),zlib=True)
    vnc.long_name = 'Latitude at vertices points'
    vnc.units = 'degrees north'
    vnc[:,:] = grdROMS.lat_vert[::decimate,::decimate]

    vnc = f1.createVariable('h', 'd', ('eta_rho','xi_rho',),zlib=True)
    vnc.long_name = 'Final bathymetry at RHO points'
    vnc.units = 'meter'
    vnc.field = "bath, scalar"
    vnc[:,:] = grdROMS.depth[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc = f1.createVariable('f', 'd', ('eta_rho','xi_rho',),zlib=True)
    vnc.long_name = 'Coriolis parameter at RHO points'
    vnc.units = 'second-1'
    vnc.field = "Coriolis, scalar"
    vnc[:,:] = grdROMS.f[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc = f1.createVariable('pm', 'd', ('eta_rho','xi_rho',),zlib=True)
    vnc.long_name = 'curvilinear coordinate metric in XI'
    vnc.units = 'meter-1'
    vnc.field = "pm, scalar"
    vnc[:,:] = grdROMS.pm[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc = f1.createVariable('pn', 'd', ('eta_rho','xi_rho',),zlib=True)
    vnc.long_name = 'curvilinear coordinate metric in ETA'
    vnc.units = 'meter-1'
    vnc.field = "pn, scalar"
    vnc[:,:] = grdROMS.pn[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc=f1.createVariable('angle','d',('eta_rho','xi_rho'),zlib=True)
    vnc.long_name = "angle between xi axis and east"
    vnc.units = "radian" 
    vnc[:,:]=grdROMS.angle[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc=f1.createVariable('mask_rho','d',('eta_rho', 'xi_rho'),zlib=True)
    vnc.long_name = "mask on RHO-points"
    vnc.option_0 = "land" 
    vnc.option_1 = "water"
    vnc.FillValue = 1.0
    vnc[:,:]=grdROMS.mask_rho[startindex:endindex:decimate,startindex:endindex:decimate]

    vnc=f1.createVariable('mask_u','d',('eta_u', 'xi_u'),zlib=True)
    vnc.long_name = "mask on U-points"
    vnc.option_0 = "land" 
    vnc.option_1 = "water"
    vnc.FillValue = 1.0
    vnc[:,:]=grdROMS.mask_u[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc=f1.createVariable('mask_v','d',('eta_v', 'xi_v'),zlib=True)
    vnc.long_name = "mask on V-points"
    vnc.option_0 = "land" 
    vnc.option_1 = "water"
    vnc.FillValue = 1.0
    vnc[:,:]=grdROMS.mask_v[startindex:endindex:decimate,startindex:endindex:decimate]
    
    vnc=f1.createVariable('mask_psi','d',('eta_psi', 'xi_psi'),zlib=True)
    vnc.long_name = "mask on PSI-points"
    vnc.option_0 = "land" 
    vnc.option_1 = "water"
    vnc.FillValue = 1.0
    vnc[:,:]=grdROMS.mask_psi[startindex:endindex:decimate,startindex:endindex:decimate]

    vnc = f1.createVariable('s_rho', 'd', ('s_rho',), zlib=True)
    vnc.long_name = "S-coordinate at RHO-points"
    vnc.valid_min = -1.
    vnc.valid_max = 0.
    vnc.field = "s_rho, scalar"
    vnc[:] = grdROMS.s_rho[:]

    vnc = f1.createVariable('s_w', 'd', ('s_w',), zlib=True)
    vnc.long_name = "S-coordinate at W-points"
    vnc.valid_min = -1.
    vnc.valid_max = 0.
    vnc.field = "s_w, scalar"
    vnc[:] = grdROMS.s_w[:]

    vnc = f1.createVariable('Cs_r', 'd', ('s_rho',), zlib=True)
    vnc.long_name = "S-coordinate stretching curves at RHO-points"
    vnc.valid_min = -1.
    vnc.valid_max = 0.
    vnc.field = "s_rho, scalar"
    vnc[:] = grdROMS.Cs_rho[:]

    vnc = f1.createVariable('Cs_w', 'd', ('s_w',), zlib=True)
    vnc.long_name = "S-coordinate stretching curves at W-points"
    vnc.valid_min = -1.
    vnc.valid_max = 0.
    vnc.field = "s_w, scalar"
    vnc[:] = grdROMS.Cs_w[:]

    vnc = f1.createVariable('hc', 'd')
    vnc.long_name = "S-coordinate parameter, critical depth";
    vnc.units = "meter"
    vnc[:] = grdROMS.hc

    vnc = f1.createVariable('Tcline', 'd')
    vnc.long_name = "S-coordinate surface/bottom layer width";
    vnc.units = "meter"
    vnc[:] = grdROMS.Tcline

    vnc = f1.createVariable('theta_s', 'd')
    vnc.long_name = "S-coordinate surface control parameter";
    vnc[:] = grdROMS.theta_s
    
    vnc = f1.createVariable('theta_b', 'd')
    vnc.long_name = "S-coordinate bottom control parameter";
    vnc[:] = grdROMS.theta_b


    f1.close()

    print "Creating new decimated grid file: %s"%(outfilename)  
    



