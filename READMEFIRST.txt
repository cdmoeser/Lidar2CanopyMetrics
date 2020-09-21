%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DESCRIPTION
%   These scripts calculates forest canopy structure metrics and density from high 
%   resolution LiDAR point cloud data of forest canopy. 

% The scripts that estimate canopy structure metrics relates to research published as
%
%     D. Moeser, F. Morsdorf, T. Jonas; Novel forest structure metrics 
%     from airborne LiDAR data for improved snow interception estimation;
%     2015; Agricultural and Forest Meteorology;
%     doi: 10.1016/j.agrformet.2015.04.013
%
%     Note: There are three versions of this script:
%
%     Version 1.1 [Lidar2CanopyMetrics] 
%	  creates a binarized canopy height raster for each evaluation point separately. 
%
%     To the contrary, version 1.2 [Lidar2CanopyMetricsAlt] 
%     creates a binarized canopy height raster encompassing the evaluation
%     domain for all evaluation points. As a consequence, version 1.2 is 
%     quicker but requires more memory. Note however, that version 1.1 creates
%     a raster which is centered around the exact coordinate of an evalution point, 
%     whereas version 1.2 uses
%     the nearest available raster cell. Therefor, paramter output from
%     both versions of the scripts is not identical, albeit very similar.
%     Version 2.0 [Lidar2CanopyMetricsGrid] 
%     is derived from version 1.2 but gives a gridded output
%     instead of point output and allows the user to choose metric(s) of
%     interest.
%
%
%     Please look at sepcific README file for version of interest:
%     1. README_Lidar2CanopyMetrics.txt
%     2. README_Lidar2CanopyMetricsAlt.txt
%     3. README_Lidar2CanopyMetricsGrid.txt
%
%     These packages come with sample DTM and LiDAR data.  
%     Note that the LIDAR folder is in .zip format due to size.  
%     If using sample data please ensure folder is unzipped
%%
%%
% The scripts that estimate canopy density metrics relates to research published as
%
%     D. Moeser, J. Roubinek, P. Schleppi, F. Morsdorf, T. Jonas;
%     Canopy closure, LAI and radiation transfer from airborne LiDAR
%     synthetic images; 2014; Agricultural and Forest Meteorology;
%     doi: 10.1016/j.agrformet.2014.06.008
   
%     Please look at sepcific README file for version of interest:
%     1. README_Lidar2CanopyMetrics.txt

% NOTE there is another reposiroty available for an updated synthetic image creator @
%     https://github.com/Tobias-Jonas-SLF/Lidar2Hemi



  
  
  
