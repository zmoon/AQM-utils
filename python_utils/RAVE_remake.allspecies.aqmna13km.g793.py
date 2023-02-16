#!/usr/bin/env python3

# Python script to handle fire emission (RAVE) to 13-km NA domain

import xarray as xr
from netCDF4 import Dataset
import ESMF
import os
import numpy as np

DATE=os.getenv('std')
year=os.getenv('yyyy')
mn=os.getenv('mn')
dd=os.getenv('dd')
cyc=os.getenv('cyc')

print(os.getenv('ravefire'))
f=os.getenv('ravefire')
ds_togid=xr.open_dataset(f)

area=ds_togid['area']
tgt_latt = ds_togid['grid_latt']
tgt_lont = ds_togid['grid_lont']

def Store_time_by_Level(fout,varname,var,long_name,yr,mm,dd,cyc,DATE1):
    if varname=='time':
        var_out = fout.createVariable(varname, 'f4', ('Time',))
        var_out.long_name = long_name
        var_out.standard_name = long_name
        fout.variables[varname][:]=var
        var_out.units = 'hours since '+yr+'-'+mm+'-'+dd+' '+cyc+':00:00'
        var_out.calendar = 'gregorian'
        var_out.axis='T'
        var_out.time_increment='010000'
        var_out.begin_date=DATE1
        var_out.begin_time='060000'

def Store_latlon_by_Level(fout,varname,var,long_name,units,dim,fval,sfactor):
    if dim=='2D':
        var_out = fout.createVariable(varname,   'f4', ('yFRP','xFRP'))
        var_out.units=units
        var_out.long_name=long_name
        var_out.standard_name=varname
        fout.variables[varname][:]=var
        var_out.FillValue=fval
        var_out.coordinates='Latitude Longitude'

def Store_by_Level(fout,varname,long_name,units,dim,fval,sfactor):
    if dim=='3D':
        var_out = fout.createVariable(varname,   'f4', ('Time','yFRP','xFRP'))
        var_out.units=units
        var_out.long_name = long_name
        var_out.standard_name=long_name
        var_out.FillValue=fval
        var_out.coordinates='Time Latitude Longitude'

outfire24h=os.getenv('outfire24h')
fout=Dataset(outfire24h,'w')
fout.createDimension('Time',24)
fout.createDimension('xFRP',800)
fout.createDimension('yFRP',544)

setattr(fout,'PRODUCT_ALGORITHM_VERSION','Beta')
setattr(fout,'TIME_RANGE','72 hours')
setattr(fout,'RangeBeginningDate\(YYYY-MM-DD\)',year+'-'+mn+'-'+dd)
setattr(fout,'RangeBeginningTime\(UTC-hour\)',cyc)
setattr(fout,'WestBoundingCoordinate\(degree\)','151.981f')
setattr(fout,'EastBoundingCoordinate\(degree\)','332.019f')
setattr(fout,'NorthBoundingCoordinate\(degree\)','81.7184f')
setattr(fout,'SouthBoundingCoordinate\(degree\)','7.22291f')

Store_latlon_by_Level(fout,'Latitude',tgt_latt,'cell center latitude','degrees_north','2D','-9999.f','1.f')
Store_latlon_by_Level(fout,'Longitude',tgt_lont,'cell center longitude','degrees_east','2D','-9999.f','1.f')

vars_emis = ["PM2.5","CO","VOCs","NOx","BC","OC","SO2","NH3","FRP_MEAN"]

for svar in vars_emis:
     if svar=='FRP_MEAN':
         Store_by_Level(fout,'MeanFRP','Mean Fire Radiative Power','MW','3D','0.f','1.f')
         tgt_rate = ds_togid[svar].fillna(0)
         fout.variables['MeanFRP'][:,:,:] = tgt_rate
     else :
         Store_by_Level(fout,svar,svar+' Biomass Emissions','kg m-2 s-1','3D','0.f','1.f')
         tgt_rate = ds_togid[svar].fillna(0)/area/3600
         fout.variables[svar][:,:,:] = tgt_rate

fout.close()
