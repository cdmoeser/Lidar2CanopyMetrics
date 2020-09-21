Lidar2CanopyMetricsGrid

%% DESCRIPTION
% This package of scripts creates forest canopy structure metrics based
% upon Aerial LiDAR and an underlying DTM.  The metrics can be
% created at each DTM grid cell or a 'sample interval' can be set which
% samples the DTM domain at a set distance. For example if the sample
% interval was set to 5, then every five meters (in the x and y direction)
% an evaluation point would be established and the algorithm would be run.
% The data derived from the vector searching algorithm is then resampled
% (interpolated) to fit the original dtm domain and cellsize.  The data is
% output as an ARC ASCII grid.
%
%PUBLICATION
%The script relates to research published as:
%   D. Moeser, F. Morsdorf, T. Jonas; Novel forest structure metrics 
%   from airborne LiDAR data for improved snow interception estimation;
%   2015; Agricultural and Forest Meteorology;
%   doi: 10.1016/j.agrformet.2015.04.013
%
%IMPLEMENTATION
%   by David Moeser and Tobias Jonas
%   WSL Institute for Snow and Avalanche Research SLF
%   Davos, Switzerland and the United States Geological Survey - New Mexico
%   Water Science Center
%
% VERSION / LAST CHANGES
%   v2.0 / 2017-08-22
%     Note: There are three versions of this script. Version 1.1 creates a 
%     binarized canopy height raster for each evaluation point separately. 
%     To the contrary, version 1.2 creates a binarized canopy height raster 
%     encompassing the evaluation domain for all evaluation points. As a 
%     consequence, version 1.2 is quicker but requires more memory. Note
%     however, that version 1.1 creates a raster which is centered around
%     the exact coordinate of an evalution point, whereas version 1.2 uses
%     the nearest available raster cell. Therefor, paramter output from
%     both versions of the scripts is not identical, albeit very similar.
%     Version 2.0 is derived from version 1.2 but gives a gridded output
%     instead of point output and allows the user to choose metric(s) of
%     interest.
% DATA REQUIREMENTS
%   Required Data Input
%      LiDAR data
%       The scripts are based around LiDAR data that comes in the standard
%       LAS format. The LAS format is a public file format for the
%       interchange of 3-dimensional point cloud data. The LAS format is
%       defined by the American Society for Photogrammetry and Remote
%       Sensing (ASPRS). Please refer to their webpage for more information
%       http://www.asprs.org/Committee-General/LASer-LAS-File-Format-Exchange-Activities.html.
%       All data parsed through this script as well as the demonstration
%       ataset use LAS 1.4. Compatibility with other versions of LAS
%       cannot be guaranteed. The LAS reader and affiliated files
%       utilized here are from Dr. Felix Morsdorf within the Remote Sensing
%       Laboratories of the University of Zürich. Of note, it is possible
%       to use other LAS readers as long as the output remains uniform to
%       the below example and there is access to the X, Y, Z data.
%     DTM data
%       For DTM data (i.e. surface elevation data without vegetation) a DTM
%       in standard geotiff or ascii format is needed.  The DTM should have a
%       footprint which is atleast as large as the input .LAS data.  The
%       DTM should have no holes.
%     PARAMETER DEFINITIONS
%       Ensure all parameters definaitions, inlcuding file paths are
%       proerly set between lines 130 - 366.
% OUTPUT
%    ARC ASCII Files(S)
%      The program will generate one ARC ASCII grid per defined metric.
%      The output grid is based upon the exact size of the .LAS file 
%      measured as a x- linear distance form the lower left corner and a 
%      y- linear distance from the lower let corner. The output file is
%      named 'metric name_lly_llx_nrows_ncols_cellsize.txt' where 'lly' is
%      follwed by the lower left corner northing, 'llx' is follwed by the
%      lower left easting
% OPTIONAL OUTPUT
%   .MAT File
%     This script will generate an output file in Matlab format, which will 
%     be saved in the output folder as specified in the user settings. The 
%     output file will be named L2CM_[yymmddTHHMM].mat according to the 
%     current time. The output file contains the struct array L2CM with the
%     following fields;
%     > L2CM.desc: cell array with a description of each output parameter;
%         rows represent different parameters
%     > L2CM.vars: data array that contains the output parameter values;  
%         rows represent different parameters, columns represent different
%         evaluation points
%     > L2CM.xcoor: data array that contains easting coordinates; columns 
%         represent different evaluation points
%     > L2CM.xcoor: data array that contains northing coordinates; columns 
%         represent different evaluation points
%     > L2CM.file: struct array with further path/filename information
%         pointing to the datasets that were used for the analysis
%     > L2CM.set: struct array that contains all settings used for the
%         analysis
%     > L2CM.version: string record on script version used
%
% CPU REQUIREMENTS
%   Note: To query a small set of points within the demonstration dataset
%   (1km x 1km) it will take approximately 2 minutes to read, parse and
%   format the .LAS data and ~15 seconds on a  current PC to parse 100
%   points for three metrics. However, as you deal with larger LiDAR
%   datasets, you may run into memory problems. In this case, consider
%   modifying the input. Please note, if it chosen to plot the output; the
%   time, CPU and RAM requirements will be significantly higher depending
%   the number of evaluation points. There are display outputs at every 10%
%   of parsed data points in order to get an ideas of run time on your
%   machine for large datasets
%
% SETUP REQUIRED TO RUN THIS SCRIPT
%   Matlab base version 7.0 or higher. As far as we are aware of no 
%   additional toolbox is needed. Add folder Lidar2CanopyMetrics 
%   including all subfolders to the Matlab path (to allow access to 
%   auxiliarys functions called from this script).
%
%  IMPORTANT NOTE       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     canopy structure metrics will be created and analyzed over the entire
%     domain. The domain is defined as the maximum extent of the current 
%     LiDAR data (.las file).  The DTM is then cut based upon the maximum 
%     extent of the curent LiDAR file. Structure metrics are calculated 
%     based upon the cut dtm at 'x' intervals. Therefore values from the 
%     canopy structure anaylsis at the DTM boundary(s) [dtm edge + 'eval_offs'] 
%     are given a value of NaN. Please consider this when inputting LiDAR 
%     data and plan accordingly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%