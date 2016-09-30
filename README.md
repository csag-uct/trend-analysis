# trend-analysis
Code related to calculating time series trends

This repo will potentially contain both Python and R versions of the algorithms but both implementations should run the same series of tests against the test data.

## Data formats

Input and output data (gridded and station) is NetCDF4 (do we want to support NetCDF3 or NetCDF4 classic?) and CF 1.6 (http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html) compliant.  Station data in NetCDF is supported under CF1.6 descrete sampling geometries (http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html#discrete-sampling-geometries) and orthogonal multi-dimensional array representation is used.

##
In the meantime, the repo contains a set of (crudely cobbled together) functions that can be used for trend analysis and derivation of trend significance level and confidence interval of individual time series, as described below:

v1.0
created by P.Wolski
September 2016

file: functions.py

functions:
test_autocorr()
trend_CI()
get_TheilSen()

 
known bugs and issues:
-lowess in trend_CI() crashes if time series starts with NAs
-there is no implementation of significance level for Durbin-Watson test for autocorrelation
-test_autocorr() and get_TheilSen() missing structured docstrings
-example datasets and wrapper functions would be nice to have
 
