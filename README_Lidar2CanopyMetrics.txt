Lidar2CanopyMetrics

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% DESCRIPTION
  %   This scripts calculates forest canopy structure metrics from high 
  %   resolution LiDAR point cloud data of forest canopy. The script
  %   relates to research published as
  %
  %     D. Moeser, F. Morsdorf, T. Jonas; Novel forest structure metrics 
  %     from airborne LiDAR data for improved snow interception estimation;
  %     2015; Agricultural and Forest Meteorology;
  %     doi: 10.1016/j.agrformet.2015.04.013
  %
  %   The script comes along with a demonstration dataset which is provided
  %   for the purpose of testing the script's functionality. Adaptations
  %   might become necessary as you apply this script to other datasets.
  %   Comments in the script may guide you through possible adaptations.
  % 
  %   Note further that this script has been developed for aerial LiDAR 
  %   data flown with a point density of approx 36pt/m2. Applications to
  %   LiDAR data from terrestrial scanning and/or with significantly
  %   different point densities may require other methodology.
  %
  % SETUP REQUIRED TO RUN THIS SCRIPT
  %   Matlab base version 7.0 or higher. As far as we are aware of no 
  %   additional toolbox is needed. Add folder Lidar2CanopyMetrics 
  %   including all subfolders to the Matlab path (to allow access to 
  %   auxiliarys functions called from this script).
  % 
  % IMPLEMENTATION
  %   by David Moeser and Tobias Jonas
  %   WSL Institute for Snow and Avalanche Research SLF
  %   Davos, Switzerland 
  %
  % VERSION / LAST CHANGES
  %   v1.2 / 2016-05-30
  %   Note: There are two versions of this script. Version 1.1 creates a 
  %   binarized canopy height raster for each evaluation point separately. 
  %   To the contrary, version 1.2 creates a binarized canopy height raster 
  %   encompassing the evaluation domain for all evaluation points. As a 
  %   consequence, version 1.2 is quicker but requires more memory. Note
  %   however, that version 1.1 creates a raster which is centered around
  %   the exact coordinate of an evalution point, whereas version 1.2 uses
  %   the nearest available raster cell. Therefore, parameter output from
  %   both versions of the scripts is not identical, albeit very similar.
  %
  % DATA REQUIREMENTS
  %   LiDAR data
  %     Our scripts are based around LiDAR data that comes is the standard
  %     LAS format. The LAS format is a public file format for the 
  %     interchange of 3-dimensional point cloud data. The LAS format is
  %     defined by the American Society for Photogrammetry and Remote
  %     Sensing (ASPRS). Please refer to their webpage for more information
  %     http://www.asprs.org/Committee-General/LASer-LAS-File-Format-Exchange-Activities.html. 
  %     All data parsed through this script as well as the demonstration
  %     dataset use LAS 1.4. Compatibility with other versions of LAS
  %     cannot be guaranteed. The LAS reader and affiliated files
  %     utilized here are from Dr. Felix Morsdorf within the Remote Sensing
  %     Laboratories of the University of Zürich. Of note, it is possible
  %     to use other LAS readers as long as the output remains uniform to
  %     the below example and there is access to the X, Y, Z data.  
  %   DTM data
  %     For DTM data (i.e. surface elevation data without vegetation) we
  %     use a list of point data. This list must consist of 3 columns for
  %     easting, northing, elevation and be saved in a text format.
  %     Refer to auxiliary function hdrload.m for more information.
  %   Evaluation points
  %     Evaluation points are coordinates for which forest canopy structure
  %     metrics are to be calculated. These points are defined from a list
  %     of point data. This list must consist of 2 columns for easting and
  %     northing and be saved in a text format.
  %   Note: all three above datasets must feature data in the same
  %     coordinate system (true coordinates). The demonstration dataset
  %     features data in the SwissGrid coordinate system CH1903/LV03.
  %
  % CPU REQUIREMENTS
  %   Note: To query a small set of points within the demonstration dataset
  %   it will take less than a minute on a current PC. However, as you deal
  %   with larger LiDAR datasets, you may run into memory problems. In this
  %   case, consider modifying the code by reading and preprocessing only 
  %   relevant data, and not the entire LiDAR / DTM files.  
  % 
  % OUTPUT
  %   This script will generate an output file in Matlab format, which will 
  %   be saved in the output folder as specified in the user settings. The 
  %   output file will be named L2CM_[yymmddTHHMM].mat according to the 
  %   current time. The output file contains the struct array L2CM with the
  %   following fields;
  %   > L2CM.desc: cell array with a description of each output parameter;
  %       rows represent different parameters
  %   > L2CM.vars: data array that contains the output parameter values;  
  %       rows represent different parameters, columns represent different
  %       evaluation points
  %   > L2CM.xcoor: data array that contains easting coordinates; columns 
  %       represent different evaluation points
  %   > L2CM.xcoor: data array that contains northing coordinates; columns 
  %       represent different evaluation points
  %   > L2CM.file: struct array with further path/filename information
  %       pointing to the datasets that were used for the analysis
  %   > L2CM.set: struct array that contains all settings used for the
  %       analysis
  %   > L2CM.version: string record on script version used
  %
  % USAGE
  %   adapt user settings as necessary
  %   run Lidar2CanopyMetrics from the command line of Matlab
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%