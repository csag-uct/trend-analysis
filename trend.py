#
# calculates trend, slope or pvalue for gridded data
#
# P.Wolski
# December 2016
#
#  
# use: trend.py input.nc output.nc linear|sen|quantreg slope|pval
#
# returns output.nc of the same lat lon dimensions as input.nc
# input.nc has to have one 3D variable with dimensions: lat,lon,time (not necessarily in this sequence)
# pval for Sen slope is calculated as pval for Mann-Kendall test
# for linear and quantreg - from analytical expressions
# 


import numpy as np
from netCDF4 import Dataset
import sys, os
import statsmodels.api as sm
from statsmodels.regression.linear_model import OLS
from statsmodels.regression.quantile_regression import QuantReg
from scipy.stats import mstats, kendalltau

#to remove duplicate warnings
import warnings
warnings.filterwarnings('ignore')


infile=sys.argv[1]
outfile=sys.argv[2]
trend=sys.argv[3]
what=sys.argv[4]

 
#



def get_TheilSen(_y, what="slope"):
    if not np.ma.is_masked(_y):
        if what=="slope":
            return mstats.theilslopes(np.ma.masked_invalid(_y))[0]
        else:
            _x=np.arange(len(_y))
            return kendalltau(_x, _y, nan_policy='omit')[1]
    return np.nan
        
def get_linear(_y, what="slope"):
    # receives 
    # data need to be regularly spaced
    # returns slope in units of _y per unit of _x, intercept in units of _y or pvalue
    # pvalue is analytical, perhaps one day I will implement bootstrap 
    #
    # need to add constant for OLS does not include intercept by default 
    if not np.ma.is_masked(_y):
        _x = sm.add_constant(np.arange(len(_y)))
        res=sm.OLS(_y, _x, missing='drop').fit()
        if what=="slope":
            return res.params[1]
        elif what=="pval":
            return res.pvalues[1]
        elif what=="intercept":
            return res.params[0]
    else:
        return np.nan

def get_quantreg(_y, what="slope", q=0.5):
    if not np.ma.is_masked(_y):    
        _x = sm.add_constant(np.arange(len(_y)))
        res=QuantReg(_y, _x).fit(q=0.5)
        if what=="slope":
            return res.params[1]
        elif what=="pval":
            return res.pvalues[1]
        elif what=="intercept":
            return res.params[0]
    else:
        return np.nan



# this reads a netCDF file with three dimensions: lon, lat and time, containing one variable
# trends are calculated on the entire scene, save the areas that have NaNs

ncdata=Dataset(infile, "r", format='NETCDF4')
#finding variable name
found=False
for v in ncdata.variables.keys():
    vdata=ncdata.variables[v]
    if v=="lat" or v=="latitude":
        latname=v
        lats=vdata[:]
        nlat=len(lats)
    if v=="lon" or v=="longitude":
        lonname=v
        lons=vdata[:]
        nlon=len(lons)
    if v=="time" or v=="date":
        timename=v
        time=vdata[:]
        nts=len(time)
    if vdata.ndim==3:
        vname=v
        data0=ncdata.variables[v][:]
        found=True
        print vname
if found==False:
    print "no appropriate variable found in the source file"
    sys.exit()

#reorder the input array, so that it is time, latitude, longitude
temp=ncdata.variables[v]
latindx=temp.dimensions.index(latname)
lonindx=temp.dimensions.index(lonname)
timeindx=temp.dimensions.index(timename)
data=np.moveaxis(data0,(timeindx,latindx,lonindx), (0,1,2))

if trend=="linear":
    res=np.ma.apply_along_axis(get_linear, 0, data, what=what)
elif trend=="TheilSen":
    res_sen=np.ma.apply_along_axis(get_TheilSen, 0, data, what=what)
elif trend=="quantile":
    res_quant=np.ma.apply_along_axis(get_quantreg, 0, data, what=what)
else:
    print "unknown trend type: "+ trend
    print "exiting..."
    sys.exit()

ncdataset=Dataset(outfile, "w", format='NETCDF3_CLASSIC')
ncdataset.createDimension('lon', nlon)
ncdataset.createDimension('lat', nlat)
longitudes = ncdataset.createVariable('lon',np.float32, ('lon',))
latitudes  = ncdataset.createVariable('lat',np.float32, ('lat',))
latitudes.long_name = "latitude" ;
latitudes.units = "degrees_north" ;
longitudes.long_name = "longitude" ;
longitudes.units = "degrees_east" ;
longitudes[:]=lons
latitudes[:]=lats

v=ncdataset.createVariable(what+"_"+trend,np.float32,('lat', 'lon',))
v[:,:]=res
if what=="slope":
    v.units="variable unit per time step"
if what=="pval":
    v.units="-"

ncdataset.close()


