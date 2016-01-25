# trend-analysis
Code related to calculating time series trends

This repo will potentially contain both Python and R versions of the algorithms but both implementations should run the same series of tests against the test data.

## Data formats

Input and output data (gridded and station) is NetCDF4 (do we want to support NetCDF3 or NetCDF4 classic?) and CF 1.6 (http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html) compliant.  Station data in NetCDF is supported under CF1.6 descrete sampling geometries (http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html#discrete-sampling-geometries) and orthogonal multi-dimensional array representation is used.
