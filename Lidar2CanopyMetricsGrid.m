function Lidar2CanopyMetricsGrid
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    %%%%%     %%%  BEGIN PARAMTER DEFINITION  %%%   %%%%%     %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    %%%%%     identify file path and file names     %%%%%     %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basefolder          = 'C:\projects\forest_snow\tilting_study\jemez\Lidar2CanopyMetricsGrid';
% path to folder which contains this script as well as all subfolders
% to auxiliary functions / data required to run this script

dtm_name         = 'DTM_data\JRB_dtm_clip.tif';
% relative path\file to data file that contains the geotiff or ascii data.
% These data must be in standard format.  For example, a standard ascii
% file should be space delimited, the first six rows contain header
% information and row 7: end contains data. The absolute path to the DTM
% data file is derived as [basefolder '\' dtm_name]

DSM                 = 'DSM_Data\JRB_points_clip.las'
% relative path\file to data file that contains DSM/nDSM data. The
% absolute path to the DSM data file is derived as [basefolder '\' DSM]

output_folder         = 'Output_Data';
% relative path to folder into which output results will be saved. The
% absolute path to the output folder is derived as
% [basefolder '\' output_folder]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    %%%%%   %%%     identify canopy metrics     %%%   %%%%%     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
METRICS         = {1; 2; 6; 10; 12; 40; 44; 48; 52; 54}
% data metrics that will be ouput in an ascii grid format.
% Please use the below list of 60 metrics and copy and paste the exact
% name(s) from the list in either the following format:
% --> METRICS = {'mean 1-sided distance to canopy'; 'total open area'; 'max tree height within local evaluation domain'}
% or  METRICS = {1; 40; 60} where the number respondes to the metrics in
% the below list

%%%% MEAN 1-SIDED DISTANCE METRICS %%%%% MEAN 1-SIDED DISTANCE METRICS %%%%
%  (1)  'mean 1-sided distance to canopy';                                 
%  (2)  'mean 1-sided distance to canopy in 90° sector to the North';   
%  (3)  'mean 1-sided distance to canopy in 90° sector to the North-East'; 
%  (4)  'mean 1-sided distance to canopy in 90° sector to the East';    
%  (5)  'mean 1-sided distance to canopy in 90° sector to the South-East'
%  (6)  'mean 1-sided distance to canopy in 90° sector to the South';   
%  (7)  'mean 1-sided distance to canopy in 90° sector to the South-West';
%  (8)  'mean 1-sided distance to canopy in 90° sector to the West';    
%  (9)  'mean 1-sided distance to canopy in 90° sector to the North-West';
%  (10) 'mean 1-sided distance to canopy in 180° sector to the North';  
%  (11) 'mean 1-sided distance to canopy in 180° sector to the East';  
%  (12) 'mean 1-sided distance to canopy in 180° sector to the South';
%  (13) 'mean 1-sided distance to canopy in 180° sector to the West'; 

%%%%% MAX 1-SIDED DISTANCE METRICS %%%%% MAX 1-SIDED DISTANCE METRICS %%%%%
%  (14) 'max 1-sided distance to canopy';       
%  (15) 'max 1-sided distance to canopy in 90° sector to the North';    
%  (16) 'max 1-sided distance to canopy in 90° sector to the North-East';
%  (17) 'max 1-sided distance to canopy in 90° sector to the East';     
%  (18) 'max 1-sided distance to canopy in 90° sector to the South-East'; 
%  (19) 'max 1-sided distance to canopy in 90° sector to the South';    
%  (20) 'max 1-sided distance to canopy in 90° sector to the South-West';
%  (21) 'max 1-sided distance to canopy in 90° sector to the West';     
%  (22) 'max 1-sided distance to canopy in 90° sector to the North-West'; 
%  (23) 'max 1-sided distance to canopy in 180° sector to the North';   
%  (24) 'max 1-sided distance to canopy in 180° sector to the East';    
%  (25) 'max 1-sided distance to canopy in 180° sector to the South';   
%  (26) 'max 1-sided distance to canopy in 180° sector to the West'; 

%%%%% MIN 1-SIDED DISTANCE METRICS %%%%% MIN 1-SIDED DISTANCE METRICS %%%%%
%  (27) 'min 1-sided distance to canopy';     
%  (28) 'min 1-sided distance to canopy in 90° sector to the North';    
%  (29) 'min 1-sided distance to canopy in 90° sector to the North-East';
%  (30) 'min 1-sided distance to canopy in 90° sector to the East';     
%  (31) 'min 1-sided distance to canopy in 90° sector to the South-East'; 
%  (32) 'min 1-sided distance to canopy in 90° sector to the South';    
%  (33) 'min 1-sided distance to canopy in 90° sector to the South-West'; 
%  (34) 'min 1-sided distance to canopy in 90° sector to the West';     
%  (35) 'min 1-sided distance to canopy in 90° sector to the North-West';
%  (36) 'min 1-sided distance to canopy in 180° sector to the North';   
%  (37) 'min 1-sided distance to canopy in 180° sector to the East'; 
%  (38) 'min 1-sided distance to canopy in 180° sector to the South';   
%  (39) 'min 1-sided distance to canopy in 180° sector to the West';  

%%%%% TOTAL OPEN AREA METRICS %%%%% TOTAL OPEN AREA METRICS %%%%% 
%  (40) 'total open area';    
%  (41) 'total open area of any 45 degree sector';                      
%  (42) 'total open area of any 90 degree sector'; 
%  (43) 'total open area area of any 180 degree sector';                
%  (44) 'total open area of 90° sector to the North';  
%  (45) 'total open area of 90° sector to the North-East';              
%  (46) 'total open area of 90° sector to the East'; 
%  (47) 'total open area of 90° sector to the South-East';              
%  (48) 'total open area of 90° sector to the South'; 
%  (49) 'total open area of 90° sector to the South-West';              
%  (50) 'total open area of 90° sector to the West'; 
%  (51) 'total open area of 90° sector to the North-West';              
%  (52) 'total open area of 180° sector to the North'; 
%  (53) 'total open area of 180° sector to the East';                   
%  (54) 'total open area of 180° sector to the South'; 
%  (55) 'total open area of 180° sector to the West';                   

%%%%% 2-SIDED DISTANCE METRICS %%%%% 2-SIDED DISTANCE METRICS %%%%%
%  (56) 'mean 2-sided distance to canopy'; 
%  (57) 'max 2-sided distance to canopy';                               
%  (58) 'min 2-sided distance to canopy';                                                                                                          

%%%%% TREE HEIGHT METRICS %%%%% TREE HEIGHT METRICS %%%%% 
%(59) 'mean tree height within local evaluation domain';
%(60) 'max tree height within local evaluation domain';               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    %%%%%   vector searching algorithm parameters     %%%%%    %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ev_offs         = 100;
% this parameter defines the evaluation domain around a point for which
% forest canopy structure metrics are to be calculated. Only LiDAR points
% within this evaluation domain will be considered when calculating the
% metrics. By default, the evaluation domain is defined as a rectangular
% box around a point with a side length of 2 ev_offs. A circular domain
% is available on request (see user setting circ_domain below). ev_offs
% is to be given in the same unit as the LiDAR data (m in case of the
% demonstration dataset). Small evaluation domains will allow for quick
% calculation, but limit the representation of extended forest gaps.
% Any extension of gaps beyond the circumference of the evaluation domain
% will be neglected.

hg_cutoff       = 1.25;
% this parameter defines a height of canopy elements above the forest
% floor underneath which all LiDAR points are neglected. A height cutoff
% is necessary to a) separate forest canopy from understory, and b) to
% avoid unrealiable output due to insufficient point cloud densities as
% you approach the forest floor (assuming LiDAR data to be taken from
% above the canopy). Suggested range for data similar to the
% demonstration dataset: 1-3.5 m.

gr_size         = 0.25;
% this parameter defines the cell size of the grid to which the LiDAR
% point cloud data will be converted. More precisly, first a tree height
% model is calculated from LiDAR and DTM data. This tree height model is
% then projected onto a rectangular grid. The cell size is to be given
% in the same unit as the LiDAR data (m in case of the demonstration
% dataset). Suggested range for data similar to the demonstration
% dataset: 0.25-1m.

negl_dist       = 0.75;
% this parameter defines a minimum distance between the above evaluation
% domain and a point for which forest canopy structure metrics are to be
% calculated. Once the tree height model has been coverted to grid, all
% grid cells that fall short of the minimum distance requirement will be
% neglected. This additional constraint has been found useful to avoid
% overinterpretation of very close canopy elements (cf. section 2.3.4
% of Moeser et al., 2015). The distance is to be given in the same unit
% as the LiDAR data (m in case of the demonstration dataset).Suggested
% range for data similar to the demonstration dataset: 0.75-1.25m.

smooth_thr      = 0.9;
% this parameter defines a threshold as part of a smoothing function
% applied to the rasterized tree height model. First, the tree height
% model is smoothed by way of a moving 3 × 3 neighborhood cell mean
% filter. The resulting grid is then re-binarized using the threshold,
% i.e. all values above or equal to the threshold are assigned to
% represent a canopy element whereas all other values are are dismissed.
% The smoothing was introduced to reduce the amount of noise effects
% from singular LiDAR returns within the lower canopy (cf. section 2.3.4
% of Moeser et al., 2015). The threshold is to be given as number within
% [0 1]. Suggested range for data similar to the demonstration dataset:
% 0.6-1.

angular_res     = 'high';
% this parameter allows to select the angular resolution which is to be
% used when calculating forest canopy structure metrics (cf. section
% 2.3.3 of Moeser et al., 2015). Chose from {'high','low'} to evaluate
% azimuth angles at a resolution of 2 and 5 degree respectively.

ev_tree_height  = 5;
% this parameter defines the evaluation domain around a point for which
% the additional metric "mean tree height" is to be calculated. Only
% LiDAR points within this evaluation domain will be considered when
% calculating this additional metric. By default, the evaluation domain
% is defined as a rectangular box around a point with a side length of 2
% ev_tree_height. A circular domain is available on request (see user
% setting circ_domain below). ev_tree_height is to be given in the same
% unit as the LiDAR data (m in case of the demonstration dataset).

circ_domain     = 1;
% set to 1 if you wish to use a circular evaluation domain, set to 0 for
% a rectangular evaluation domain (as in Moeser et al., 2015)

doplot          = 0;
% set to 1 if you wish to display a plot for the initiual domain with the
% polygons derived from the vector seraching algorithm from all analysis
% points plotted on top, set to 0 otherwise. The plot will show an analysis
% of the gap space around each evaluation point. Please note, if it chosen 
% to plot the output, the time, CPU and RAM requirements will be 
% significantly higher depending the number of evaluation points. 

sversion        = 'v2.0  / 2017-08-25';
% string, will be save in the results file for reference purposes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    %%%%%     define geotiff data properties     %%%%%     %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nodata               = NaN;
%define nodata. This replaces any values less than 0 to this value.  It is
%reocmended this be set to NaN
cellsize             = 1;
%define native cellsize of geotiff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    %%%%%     define sample interval    %%%%%     %%%%%%     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sample_interval     = 100;
% define every 'x' meters a syntheic image is to be created and analyzed.
% This should be larger than the 'cellsize'  If not this values defaults
% to 'cellsize'
interpolation_method = 'linear'; 
% define interpolation method to return to native geotiff cell size.
% nearest neighbor or 'nearest', 'linear interpolation or linear', or
% natural neighbor 'natural' can be used  if blank ('') this method defaults
% to 'linear'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    %%%%%     %% optional data output %%   %%%%%     %%%%%%     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meta_raw_data_save =  1;
% set to 1 if you want to save the meta data and raw data in a matlab
% structure format (saved in base folder) where:
%     L2CM.desc                   (parameter description)
%     L2CM.vars                   (parameter values)
%     L2CM.xcoor                  (easting of evaluation points)
%     L2CM.ycoor                  (northing of evaluation points)
%     L2CM.file.path              (files used)
%     L2CM.file.LiDAR             (DSM file used)
%     L2CM.file.DTM               (DTM file used)
%     L2CM.set.eval_max_dist      (settings used - max evaluation distance)
%     L2CM.set.eval_min_dist      (settings used - min evaluation distance)
%     L2CM.set.eval_loc_dist      (settings used - evaluation domain for tree heigh metrics)
%     L2CM.set.eval_shape         (settings used - evaluation domain shape -> circular or rectangular)
%     L2CM.set.height_cutoff      (settings used - height of canopy elements above the forest floor underneath which all LiDAR points are neglected)
%     L2CM.set.raster_cellsize    (settings used - cell size of the converted LiDAR point cloud)
%     L2CM.set.raster_smoothing   (settings used - smoothing function value)
%     L2CM.set.angular_resolution (settings used - angualr resolution of the vector search --> 2 or 5 degrees)
%     L2CM.version                (program version)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    %%%%%     %%%  END PARAMTER DEFINITION    %%%   %%%%%     %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    %%%%%     %%%%%     PREPARATION STEPS     %%%%%   %%%%%     %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic % start timer for display messages
  %%% add path to aux functions 
  addpath([basefolder '\Aux_Functions']);
%% read header information form LAS data to collect data extents
out          = readlas_hdr([basefolder '\' DSM] );
max_easting  = out.MaxX;
min_easting  = out.MinX;
max_northing = out.MaxY;
min_northing = out.MinY;
easting      = [min_easting max_easting];   % put data into correct format for geoimread2
northing     = [min_northing max_northing];
%%  read geotiff or ascii data based on the extent of the .LAS data and convert nodata values to 'nodata' values  and create
% ascii data must be the STANDARD format (space delineated, first six rows
% contain header information.
% (a) vector list of x, y, data based upon actual geotiff spacing (cellsize)
% (b) vector list of x, y, data based upon sample interval (sample_interval) 

if contains(dtm_name, '.tif') ==1
    dtm_data                 = geoimread2([basefolder '\' dtm_name], easting, northing);   % [dtm_data, x, y, I]  = geoimread2([geotiff_path dtm_name], easting, northing); also gives controid x (x), ceontroid y values (y), and dtm_data projection information (I)
elseif contains(dtm_name, '.tif') ~=1
    dtm_data = dlmread([basefolder '\' dtm_name],' ', 6,0);
end
dtm_data(dtm_data < 0 )   = nodata;                                                         % change NaN data format to nodata value;    
    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     START 'A'     %%%%%     START 'A'     %%%%%     START 'A'     %%%%%     START 'A'     %%%%%     START 'A'     %%%%%     START 'A'     %%%%%     START 'A'     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start on the lower left corner where y and x is minimum
% create x y z vectors of data based upon the center of each grid cell in the dtm_data read-in
offset_actual = cellsize/2;            % create offset based upon cell size
ncols         = length(dtm_data(1,:));  % get number of columns (x or easting) based upon dtm_data
nrows         = length(dtm_data(:,1));  % get number of rows (y or northing) based upon dtm_data

xv            = (min_easting:cellsize:min_easting+((ncols-1)*cellsize))+offset_actual;       % get center x values (easting values)
dtm.x          = double(repmat(xv,[1 nrows])');                                              % reapeat values based upon number of rows (number of northing values) and create vector

yv            = ((min_northing :cellsize:min_northing+((nrows-1)*cellsize))+offset_actual); % get center y values (northing values)
dtm.y          = repmat(yv,[ncols 1]);                                                      % reapeat values based upon number of rows (number of easting values) and create vector
dtm.y          = double(dtm.y(:));

dtm.z          = fliplr(dtm_data');                                                          % get z values based upon the dtm_data
dtm.z          = double(dtm.z(:));                                                          % create vector   % this will return to original matrix  -->  =  Z = flipud(reshape(dtm.z, [ncols,nrows])');

clear xv  yv
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     START 'B'     %%%%%     START 'B'     %%%%%     START 'B'     %%%%%     START 'B'     %%%%%     START 'B'     %%%%%     START 'B'     %%%%%     START 'B'     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   !!!!!!!!   -->> If the actual number of nrows or ncols are not divisble by the 'sample_interval' then the number is rounded down to the nearest whole number
%   !!!!!!!!   -->> This ensures that no sampling takes place outside of the bounds of the LiDAR data
sample_nrows  = floor(nrows / sample_interval);
sample_ncols  = floor(ncols / sample_interval);

offset_sample = sample_interval/2;            % create offset based upon 'sample_interval'

xv            = (min_easting:sample_interval:min_easting+((sample_ncols-1)*sample_interval))+offset_sample;       % get center x values (easting values)
xvec          = repmat(xv,[1 sample_nrows])';                                                                     % reapeat values based upon number of rows (number of northing values) and create vector

yv            = ((min_northing :sample_interval:min_northing+((sample_nrows-1)*sample_interval))+offset_sample);  % get center y values (northing values)
yvec          = repmat(yv,[sample_ncols 1]);                                                                      % reapeat values based upon number of columns (number of easting values)
yvec          = yvec(:);                                                                                          % create vector

clear xv yv

% --> get coordinate of points for which synthethic hemispherical images are to be calculated along with aspect and gradient data based upon sample spacing % 'sample_interval'
ev_coors = double([xvec yvec]);

clear  xvec yvec xv_relative xvec_relative yv_relative yvec_relative yvec_relative2 newvec step Znew grad asp out max_easting min_easting max_northing min_northing
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     START 'C'     %%%%%     START 'C'     %%%%%     START 'C'     %%%%%     START 'C'     %%%%%     START 'C'     %%%%%     START 'C'     %%%%%     START 'C'     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% > read DSM data
dsm           = readlas([basefolder '\' DSM]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% > preprocess DSM/DTM data

%%% interpolate DTM to DSM coordinates (point to point)
si    = scatteredInterpolant(dtm.x,dtm.y,dtm.z);
dsm.e = si(dsm.x,dsm.y);
%%% convert DSM to nDSM if not done already
if mode(dsm.z) > 100
    dsm.z = dsm.z - dsm.e;    % height above terrain
end

%%% throw away z values outside dtm and hth (if raw.z is less than height cutoff then throw out value)
dsm.z(dsm.z<=hg_cutoff) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% > prepare binarized canopy height raster

%%% create raster based upon 'gr_size' from initial 'dtm_data'


[thr.xm,thr.ym]       = meshgrid(easting(1):gr_size:easting(2),northing(1):gr_size:northing(2)); % defines grid coordinates
nisix                 = find(~isnan(dsm.z));


si2                  = scatteredInterpolant(dsm.x(nisix),dsm.y(nisix),dsm.z(nisix));
thr.zm                = si2(thr.xm,thr.ym);

%%% binarize canopy height raster
thr.zb                = thr.zm;
thr.zb(thr.zm < hg_cutoff) = 0;
thr.zb(thr.zm >= hg_cutoff) = 1;

%%% smooth binarized canopy height raster using a moving 3 × 3 neighborhood cell mean filter
thr.zb = smooth2a(thr.zb,1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% define pointers to search directions (azimuth angles, only approximately equiangular approach)
switch angular_res
    case 'low' % search for canopy elements every 5° (72 directions)
        dir_vec = [0 1; 0.111 1; 0.222 1; 0.333 1; 0.444 1; 0.555 1; 0.666 1; 0.777 1; 0.888	1; 1 1; 1 0.888; 1	0.777; 1 0.666; 1 0.555; 1 0.444; 1	0.333; 1 0.222; 1 0.111; 1	0; 1 -0.111; 1 -0.222;...
            1 -0.333; 1 -0.444; 1 -0.555; 1 -0.666; 1 -0.777; 1 -0.888; 1 -1; 0.888 -1; 0.777 -1; 0.666 -1; 0.555 -1; 0.444 -1; 0.333 -1; 0.222 -1; 0.111 -1; 0 -1; -0.111 -1; -0.222 -1; -0.333	-1;...
            -0.444 -1; -0.555 -1; -0.666 -1; -0.777	-1; -0.888 -1; -1 -1; -1 -0.888; -1	-0.777; -1 -0.666; -1 -0.555; -1 -0.444; -1	-0.333; -1 -0.222; -1 -0.111; -1 0; -1 0.111; -1 0.222; -1 0.333;...
            -1 0.444; -1 0.555; -1 0.666; -1 0.777; -1 0.888; -1	1; -0.888 1; -0.777	1; -0.666 1; -0.555	1; -0.444 1; -0.333	1; -0.222 1; -0.111	1];
    case 'high' % search for canopy elements every 2° (192 directions)
        dir_vec = [0.0416666666666667 1;	0.0833333333333333 1;	0.125 1;	0.166666666666667 1;	0.208333333333333 1;	0.25 1;	0.291666666666667 1;	0.333333333333333 1;	0.375 1;	0.416666666666667 1;	0.458333333333333 1;	0.5 1;	0.541666666666667 1;	0.583333333333333 1;...
            0.625 1;	0.666666666666667 1;	0.708333333333333 1;	0.75 1;	0.791666666666667 1;	0.833333333333333 1;	0.875 1;	0.916666666666667 1;	0.958333333333333 1;	1 1;	1 0.958333333333333;	1 0.916666666666667;	1 0.875;	1 0.833333333333333;	1 0.791666666666667;...
            1 0.75;	1 0.708333333333333;	1 0.666666666666667;	1 0.625;	1 0.583333333333333;	1 0.541666666666667;	1 0.5;	1 0.458333333333333;	1 0.416666666666667;	1 0.375;	1 0.333333333333333;	1 0.291666666666667;	1 0.25;	1 0.208333333333333;	1 0.166666666666667;...
            1 0.125;	1 0.0833333333333334;	1 0.0416666666666666;	1 0;	1 -0.0416666666666667;	1 -0.0833333333333333;	1 -0.125;	1 -0.166666666666667;	1 -0.208333333333333;	1 -0.25;	1 -0.291666666666667;	1 -0.333333333333333;	1 -0.375;	1 -0.416666666666667;	1 -0.458333333333333;...
            1 -0.5;	1 -0.541666666666667;	1 -0.583333333333333;	1 -0.625;	1 -0.666666666666667;	1 -0.708333333333333;	1 -0.75;	1 -0.791666666666667;	1 -0.833333333333333;	1 -0.875;	1 -0.916666666666667;	1 -0.958333333333333;	1 -1;	0.958333333333333 -1;	0.916666666666667 -1;	0.875 -1;...
            0.833333333333333 -1;	0.791666666666667 -1;	0.75 -1;	0.708333333333333 -1;	0.666666666666667 -1;	0.625 -1;	0.583333333333333 -1;	0.541666666666667 -1;	0.5 -1;	0.458333333333333 -1;	0.416666666666667 -1;	0.375 -1;	0.333333333333333 -1;	0.291666666666667 -1;	0.25 -1;	0.208333333333333 -1;...
            0.166666666666667 -1;	0.125 -1;	0.0833333333333334 -1;	0.0416666666666666 -1;	0 -1;	-0.0416666666666667 -1;	-0.0833333333333333 -1;	-0.125 -1;	-0.166666666666667 -1;	-0.208333333333333 -1;	-0.25 -1;	-0.291666666666667 -1;	-0.333333333333333 -1;	-0.375 -1;	-0.416666666666667 -1;...
            -0.458333333333333 -1;	-0.5 -1;	-0.541666666666667 -1;	-0.583333333333333 -1;	-0.625 -1;	-0.666666666666667 -1;	-0.708333333333333 -1;	-0.75 -1;	-0.791666666666667 -1;	-0.833333333333333 -1;	-0.875 -1;	-0.916666666666667 -1;	-0.958333333333333 -1;	-1 -1;	-1 -0.958333333333333;...
            -1 -0.916666666666667;	-1 -0.875;	-1 -0.833333333333333;	-1 -0.791666666666667;	-1 -0.75;	-1 -0.708333333333333;	-1 -0.666666666666667;	-1 -0.625;	-1 -0.583333333333333;	-1 -0.541666666666667;	-1 -0.5;	-1 -0.458333333333333;	-1 -0.416666666666667;	-1 -0.375;	-1 -0.333333333333333;...
            -1 -0.291666666666667;	-1 -0.25;	-1 -0.208333333333333;	-1 -0.166666666666667;	-1 -0.125;	-1 -0.0833333333333334;	-1 -0.0416666666666666;	-1 0;	-1 0.0416666666666667;	-1 0.0833333333333333;	-1 0.125;	-1 0.166666666666667;	-1 0.208333333333333;	-1 0.25;	-1 0.291666666666667;	-1 0.333333333333333;...
            -1 0.375;	-1 0.416666666666667;	-1 0.458333333333333;	-1 0.5;	-1 0.541666666666667;	-1 0.583333333333333;	-1 0.625;	-1 0.666666666666667;	-1 0.708333333333333;	-1 0.75;	-1 0.791666666666667;	-1 0.833333333333333;	-1 0.875;	-1 0.916666666666667;	-1 0.958333333333333;	-1 1;	-0.958333333333333 1;	-0.916666666666667 1;...
            -0.875 1;	-0.833333333333333 1;	-0.791666666666667 1;	-0.75 1;	-0.708333333333333 1;	-0.666666666666667 1;	-0.625 1;	-0.583333333333333 1;	-0.541666666666667 1;	-0.5 1;	-0.458333333333333 1;	-0.416666666666667 1;	-0.375 1;	-0.333333333333333 1;	-0.291666666666667 1;	-0.25 1;	-0.208333333333333 1;	-0.166666666666667 1;...
            -0.125 1;	-0.0833333333333334 1;	-0.0416666666666666 1;	0 1;];
    otherwise
        error('Error, unknown choice as to angular resolution')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc

% create stauts updates
disp(['  -->  Done Reading and Preparing Point Cloud Data with an analysis time of ' num2str(toc/60) ' Minutes  <--']);   % create stauts bar updates
disp('  -->  Initiating Vector Searching Algorithm  <--');
disp(['  --> For a Total of ' num2str(length(ev_coors(:,1))) ' Points'])
disp('  -->  There Will Be Status Updates For Every 10 Percent of Analysis Points  <--');
status_percent = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]; % create stauts updates based upon percentage
%% >  loop through coordinates of evaluation points
tic
for coordidx = 1:length(ev_coors(:,1))                      % loop through coordinates of evaluation points
    
    % create stauts updates based upon percentage
    for sp_i = 1:length(status_percent)
        if coordidx == round(status_percent(sp_i)*(length(ev_coors(:,1))))
            disp (' --> THINGS ARE GOING REALLY WELL  <--')
            disp ([' --> THE PROGRAM HAS ANALYZED ' num2str(status_percent(sp_i)*100) ' PERCENT OF ALL POINTS or (' num2str(round(status_percent(sp_i)*(length(ev_coors(:,1))))) ' points)  <--'])
            disp ([' --> WITH A TOTAL TIME OF ' num2str(toc/60) ' MINUTES OR ' num2str(toc/3600) ' HOURS  <--'])
        end
    end
    %%% create canopy height raster for evaluation domain
    x_coor            = ev_coors(coordidx,1);
    y_coor            = ev_coors(coordidx,2);
    
    %%% identify position of evaluation point within canopy height raster
    nox               = length(thr.xm(1,:));  % number of columns
    noy               = length(thr.ym(:,1));  % number of rows
    xll               = max(thr.xm(1,thr.xm(1,1:end) <= x_coor));
    yll               = max(thr.ym(thr.ym(1:end,1) <= y_coor,1));
    if isempty(xll) || x_coor-xll > gr_size || isempty(yll) || y_coor-yll > gr_size
        error('Error locating position of point for which forest canopy structure metrics are to be calculated within canopy height raster, this should not happen')
    else
        xix = find(abs(thr.xm(1,:)-xll)<eps);
        yix = find(abs(thr.ym(:,1)-yll)<eps);
    end
    negl_ix = round(negl_dist/gr_size);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% > analyze gap space around evaluation point
    
    poly_x = [];
    poly_y = [];
    poly_d = [];
    
    for dirix = 1:size(dir_vec,1) % loop through directions [N, NNE, NE, ENE, E, ESE, SE, SSE, S, SSW, SW, WSW, W, WNW, NW, NNW, ETC --> based upon 'dir_vec'
        
        
        
        % if the coordinate is less than the 'ev_off' distance away
        % from the domain border then it is given NaN if greater than
        % the ev_off distnace then calculations begin
        if  (easting(1) + ev_offs) > x_coor ||   (easting(2) - ev_offs) < x_coor || (northing(1) + ev_offs) > y_coor || (northing(2) - ev_offs) < y_coor
            poly_d(dirix) = NaN;
            poly_x(dirix) = NaN;
            poly_y(dirix) = NaN;
        else
            
            %%% start at evaluation point and move forward into search direction by minimum distance
            cposix = round([xix+negl_ix*dir_vec(dirix,1),yix+negl_ix*dir_vec(dirix,2)]);
            if cposix(1) > nox || cposix(1) < 1 || cposix(2) > noy || cposix(2) < 1;
                % condition met if evaluation point is outside evaluation domain (should not happen)
                % anyway, in this case assume hitting a canopy element
                poly_d(dirix) = (sqrt((cposix(1)-xix)^2 + (cposix(2)-yix)^2))*gr_size;
                poly_x(dirix) = x_coor+(cposix(1)-xix)*gr_size;
                poly_y(dirix) = y_coor+(cposix(2)-yix)*gr_size;
            elseif thr.zb(cposix(2),cposix(1)) >= smooth_thr;
                % condition met if a tree is hit
                poly_d(dirix) = (sqrt((cposix(1)-xix)^2 + (cposix(2)-yix)^2))*gr_size;
                poly_x(dirix) = x_coor+(cposix(1)-xix)*gr_size;
                poly_y(dirix) = y_coor+(cposix(2)-yix)*gr_size;
            else
                
                %%% iteratively move forward into search direction until a tree / or the boundary of the evaluation domain is hit
                nohit        = 1;
                step_count   = negl_ix;
                while nohit;
                    step_count = step_count+1;
                    cposix     = round([xix+step_count*dir_vec(dirix,1),yix+step_count*dir_vec(dirix,2)]);
                    if cposix(1) > nox || cposix(1) < 1 || cposix(2) > noy || cposix(2) < 1;
                        % condition met if current position is outside evaluation domain
                        % in this case assume hitting a canopy element
                        poly_d(dirix) = (sqrt((cposix(1)-xix)^2 + (cposix(2)-yix)^2))*gr_size;
                        poly_x(dirix) = x_coor+(cposix(1)-xix)*gr_size;
                        poly_y(dirix) = y_coor+(cposix(2)-yix)*gr_size;
                        nohit = 0;
                    elseif thr.zb(cposix(2),cposix(1)) >= smooth_thr || isnan(thr.zb(cposix(2),cposix(1)));
                        % condition met if a tree is hit
                        poly_d(dirix) = (sqrt((cposix(1)-xix)^2 + (cposix(2)-yix)^2))*gr_size;
                        poly_x(dirix) = x_coor+(cposix(1)-xix)*gr_size;
                        poly_y(dirix) = y_coor+(cposix(2)-yix)*gr_size;
                        nohit = 0;
                    end
                end
            end
        end
    end
    
    %%% optionally restrict results to circular domain
    if circ_domain
        for dirix = 1:size(dir_vec,1)
            if poly_d(dirix) > ev_offs
                poly_x(dirix) = x_coor+(poly_x(dirix)-x_coor)*ev_offs/poly_d(dirix);
                poly_y(dirix) = y_coor+(poly_y(dirix)-y_coor)*ev_offs/poly_d(dirix);
                poly_d(dirix) = ev_offs;
            end
        end
    end
    
    %%% create a plot of the domain with the vector searching  polygon of
    %%% the gap space around the evaluation points displayed
    if coordidx == 1  % create the base figure
        if doplot
            fh = figure;
            ph = pcolor(thr.xm,thr.ym,thr.zb); hold all
            set(ph,'linestyle','none');
            plot(x_coor, y_coor,'.K');
            pg = fill(poly_x+1/2*gr_size,poly_y+1/2*gr_size,'w');
            set(pg,'EdgeColor',[0 1 1],'FaceColor','none','LineWidth',2)
            for dirix = 1:size(dir_vec,1)
                plot([x_coor,poly_x(dirix)]+1/2*gr_size,[y_coor,poly_y(dirix)]+1/2*gr_size,'g-');
            end
            set(gca,'xlim',[min(thr.xm(~isnan(thr.zb))),max(thr.xm(~isnan(thr.zb)))],'ylim',[min(thr.ym(~isnan(thr.zb))),max(thr.ym(~isnan(thr.zb)))]);
            drawnow;
        end
        
    elseif coordidx > 1  % add the polygons per evaluation point
        if doplot
        pg = fill(poly_x+1/2*gr_size,poly_y+1/2*gr_size,'w');
        set(pg,'EdgeColor',[0 1 1],'FaceColor','none','LineWidth',2)
        for dirix = 1:size(dir_vec,1)
            plot([x_coor,poly_x(dirix)]+1/2*gr_size,[y_coor,poly_y(dirix)]+1/2*gr_size,'g-');
        end
        set(gca,'xlim',[min(thr.xm(~isnan(thr.zb))),max(thr.xm(~isnan(thr.zb)))],'ylim',[min(thr.ym(~isnan(thr.zb))),max(thr.ym(~isnan(thr.zb)))]);
        drawnow;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% > calculate / compile forest canopy structure variables
    fcs_vars = [];
    for m_i = 1:length(METRICS)
        
        if strcmpi('mean 2-sided distance to canopy', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  &&  METRICS{m_i} == 56)
            switch angular_res
                case 'low'
                    fcs_vars(m_i)  = mean(poly_d(1:36)+poly_d(37:72));
                case 'high'
                    fcs_vars(m_i)  = mean(poly_d(1:96)+poly_d(97:192));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('max 2-sided distance to canopy', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}==57)
            switch angular_res
                case 'low'
                    fcs_vars(m_i)  = max(poly_d(1:36)+poly_d(37:72));
                case 'high'
                    fcs_vars(m_i)  = max(poly_d(1:96)+poly_d(97:192));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('min 2-sided distance to canopy', METRICS{m_i})==1|| ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 58)
            switch angular_res
                case 'low'
                    fcs_vars(m_i)  = min(poly_d(1:36)+poly_d(37:72));
                case 'high'
                    fcs_vars(m_i)  = min(poly_d(1:96)+poly_d(97:192));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('mean 1-sided distance to canopy', METRICS{m_i})==1|| ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 1)
            switch angular_res
                case 'low'
                    fcs_vars(m_i)  = mean(poly_d(1:72));
                case 'high'
                    fcs_vars(m_i)  = mean(poly_d(1:192));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('max 1-sided distance to canopy', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 14)
            switch angular_res
                case 'low'
                    fcs_vars(m_i)  = max(poly_d(1:72));
                case 'high'
                    fcs_vars(m_i)  = max(poly_d(1:192));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('min 1-sided distance to canopy', METRICS{m_i})==1|| ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 27)
            switch angular_res
                case 'low'
                    fcs_vars(m_i)  = min(poly_d(1:72));
                case 'high'
                    fcs_vars(m_i)  = min(poly_d(1:192));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('total open area', METRICS{m_i})==1|| ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 40)
            switch angular_res
                case 'low'
                    fcs_vars(m_i)  = polyarea(poly_x,poly_y);
                case 'high'
                    fcs_vars(m_i)  = polyarea(poly_x,poly_y);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
            %%% parameters calculated over moving angular sectors
        elseif strcmpi('total open area of any 45 degree sector', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 41)
            switch angular_res
                case 'low'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_45        = [];
                    for dirix = 1:size(dir_vec,1)
                        pa_45(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+9)],[y_coor,poly_yy(dirix:dirix+9)]);
                    end
                    fcs_vars(m_i)  = max(pa_45);
                case 'high'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_45        = [];
                    for dirix = 1:size(dir_vec,1)
                        pa_45(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+24)],[y_coor,poly_yy(dirix:dirix+24)]);
                    end
                    fcs_vars(m_i)  = max(pa_45);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('total open area of any 90 degree sector', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 42)
            switch angular_res
                case 'low'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_90        = [];
                    for dirix = 1:size(dir_vec,1)
                        pa_90(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+18)],[y_coor,poly_yy(dirix:dirix+18)]);
                    end
                    fcs_vars(m_i)  = max(pa_90);
                case 'high'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_90        = [];
                    for dirix = 1:size(dir_vec,1)
                        pa_90(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+48)],[y_coor,poly_yy(dirix:dirix+48)]);
                    end
                    fcs_vars(m_i)  = max(pa_90);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('total open area of any 180 degree sector', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 43)
            switch angular_res
                case 'low'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_180       = [];
                    for dirix = 1:size(dir_vec,1)
                        pa_180(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+36)],[y_coor,poly_yy(dirix:dirix+36)]);
                    end
                    fcs_vars(m_i)  = max(pa_180);
                case 'high'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_180       = [];
                    for dirix = 1:size(dir_vec,1)
                        pa_180(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+96)],[y_coor,poly_yy(dirix:dirix+96)]);
                    end
                    fcs_vars(m_i)  = max(pa_180);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
            %%% parameters calculated over fix 90° sectors
        elseif strcmpi('total open area of 90° sector to the North', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 44)
            switch angular_res
                case 'low'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_90        = [];
                    for dirix = 1:size(dir_vec,1)
                        pa_90(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+18)],[y_coor,poly_yy(dirix:dirix+18)]);
                    end
                    fcs_vars(m_i)  = pa_90(64);
                case 'high'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_90        = [];
                   for dirix = 1:size(dir_vec,1)
                        pa_90(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+48)],[y_coor,poly_yy(dirix:dirix+48)]);
                   end;
                    fcs_vars(m_i)  = pa_90(168);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('mean 1-sided distance to canopy in 90° sector to the North', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 2)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = mean([poly_d(64:72),poly_d(1:10)]);
                case 'high'
                    fcs_vars(m_i) = mean([poly_d(168:192),poly_d(1:24)]);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('max 1-sided distance to canopy in 90° sector to the North', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 15)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = max([poly_d(64:72),poly_d(1:10)]);
                case 'high'
                    fcs_vars(m_i) = max([poly_d(168:192),poly_d(1:24)]);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('min 1-sided distance to canopy in 90° sector to the North', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 28)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = min([poly_d(64:72),poly_d(1:10)]);
                case 'high'
                    fcs_vars(m_i) = min([poly_d(168:192),poly_d(1:24)]);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('total open area of 90° sector to the North-East', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 45)
            switch angular_res
                case 'low'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_90        = [];
                    for dirix = 1:size(dir_vec,1)
                        pa_90(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+18)],[y_coor,poly_yy(dirix:dirix+18)]);
                    end
                    fcs_vars(m_i) = pa_90(1);
                case 'high'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_90        = [];
                   for dirix = 1:size(dir_vec,1)
                        pa_90(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+48)],[y_coor,poly_yy(dirix:dirix+48)]);
                   end;
                    fcs_vars(m_i) = pa_90(1);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('mean 1-sided distance to canopy in 90° sector to the North-East', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 3)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = mean(poly_d(1:19));
                case 'high'
                    fcs_vars(m_i) = mean(poly_d(1:48));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('max 1-sided distance to canopy in 90° sector to the North-East', METRICS{m_i})==1 ||  ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 16)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = max(poly_d(1:19));
                case 'high'
                    fcs_vars(m_i) = max(poly_d(1:48));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('min 1-sided distance to canopy in 90° sector to the North-East', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 29)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = min(poly_d(1:19));
                case 'high'
                    fcs_vars(m_i) = min(poly_d(1:48));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('total open area of 90° sector to the East', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 46)
            switch angular_res
                case 'low'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_90        = [];
                    for dirix = 1:size(dir_vec,1)
                        pa_90(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+18)],[y_coor,poly_yy(dirix:dirix+18)]);
                    end
                    fcs_vars(m_i) = pa_90(10);
                case 'high'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_90        = [];
                   for dirix = 1:size(dir_vec,1)
                        pa_90(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+48)],[y_coor,poly_yy(dirix:dirix+48)]);
                   end;
                    fcs_vars(m_i) = pa_90(24);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('mean 1-sided distance to canopy in 90° sector to the East', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 4)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = mean(poly_d(10:28));
                case 'high'
                    fcs_vars(m_i) = mean(poly_d(24:72));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('max 1-sided distance to canopy in 90° sector to the East', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 17)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = max(poly_d(10:28));
                case 'high'
                    fcs_vars(m_i) = max(poly_d(24:72));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('min 1-sided distance to canopy in 90° sector to the East', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 30)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = min(poly_d(10:28));
                case 'high'
                    fcs_vars(m_i) = min(poly_d(24:72));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('total open area of 90° sector to the South-East', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 47)
            switch angular_res
                case 'low'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_90        = [];
                    for dirix = 1:size(dir_vec,1)
                        pa_90(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+18)],[y_coor,poly_yy(dirix:dirix+18)]);
                    end
                    fcs_vars(m_i) = pa_90(19);
                case 'high'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_90        = [];
                   for dirix = 1:size(dir_vec,1)
                        pa_90(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+48)],[y_coor,poly_yy(dirix:dirix+48)]);
                   end;
                    fcs_vars(m_i) = pa_90(48);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('mean 1-sided distance to canopy in 90° sector to the South-East', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 5)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = mean(poly_d(19:37));
                case 'high'
                    fcs_vars(m_i) = mean(poly_d(48:96));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('max 1-sided distance to canopy in 90° sector to the South-East', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 18)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = max(poly_d(19:37));
                case 'high'
                    fcs_vars(m_i) = max(poly_d(48:96));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('min 1-sided distance to canopy in 90° sector to the South-East', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 31)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = min(poly_d(19:37));
                case 'high'
                    fcs_vars(m_i) = min(poly_d(48:96));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('total open area of 90° sector to the South', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 48)
            switch angular_res
                case 'low'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_90        = [];
                    for dirix = 1:size(dir_vec,1)
                        pa_90(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+18)],[y_coor,poly_yy(dirix:dirix+18)]);
                    end
                    fcs_vars(m_i) = pa_90(28);
                case 'high'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_90        = [];
                   for dirix = 1:size(dir_vec,1)
                        pa_90(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+48)],[y_coor,poly_yy(dirix:dirix+48)]);
                   end;
                    fcs_vars(m_i) = pa_90(72);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('mean 1-sided distance to canopy in 90° sector to the South', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 6)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = mean(poly_d(28:46));
                case 'high'
                    fcs_vars(m_i) = mean(poly_d(72:120));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('max 1-sided distance to canopy in 90° sector to the South', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 19)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = max(poly_d(28:46));
                case 'high'
                    fcs_vars(m_i) = max(poly_d(72:120));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('min 1-sided distance to canopy in 90° sector to the South', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 32)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = min(poly_d(28:46));
                case 'high'
                    fcs_vars(m_i) = min(poly_d(72:120));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('total open area of 90° sector to the South-West', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 49)
            switch angular_res
                case 'low'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_90        = [];
                    for dirix = 1:size(dir_vec,1)
                        pa_90(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+18)],[y_coor,poly_yy(dirix:dirix+18)]);
                    end
                    fcs_vars(m_i) = pa_90(37);
                case 'high'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_90        = [];
                   for dirix = 1:size(dir_vec,1)
                        pa_90(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+48)],[y_coor,poly_yy(dirix:dirix+48)]);
                   end;
                    fcs_vars(m_i) = pa_90(96);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('mean 1-sided distance to canopy in 90° sector to the South-West', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 7)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = mean(poly_d(37:55));
                case 'high'
                    fcs_vars(m_i) = mean(poly_d(96:144));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('max 1-sided distance to canopy in 90° sector to the South-West', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 20)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = max(poly_d(37:55));
                case 'high'
                    fcs_vars(m_i) = max(poly_d(96:144));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('min 1-sided distance to canopy in 90° sector to the South-West', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 33)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = min(poly_d(37:55));
                case 'high'
                    fcs_vars(m_i) = min(poly_d(96:144));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('total open area of 90° sector to the West', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 50)
            switch angular_res
                case 'low'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_90        = [];
                    for dirix = 1:size(dir_vec,1)
                        pa_90(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+18)],[y_coor,poly_yy(dirix:dirix+18)]);
                    end
                    fcs_vars(m_i) = pa_90(46);
                case 'high'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_90        = [];
                   for dirix = 1:size(dir_vec,1)
                        pa_90(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+48)],[y_coor,poly_yy(dirix:dirix+48)]);
                   end;
                    fcs_vars(m_i) = pa_90(120);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('mean 1-sided distance to canopy in 90° sector to the West', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 8)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = mean(poly_d(46:64));
                case 'high'
                    fcs_vars(m_i) = mean(poly_d(120:168));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('max 1-sided distance to canopy in 90° sector to the West', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 21)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = max(poly_d(46:64));
                case 'high'
                    fcs_vars(m_i) = max(poly_d(120:168));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('min 1-sided distance to canopy in 90° sector to the West', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 34)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = min(poly_d(46:64));
                case 'high'
                    fcs_vars(m_i) = min(poly_d(120:168));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('total open area of 90° sector to the North-West', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 51)
            switch angular_res
                case 'low'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_90        = [];
                    for dirix = 1:size(dir_vec,1)
                        pa_90(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+18)],[y_coor,poly_yy(dirix:dirix+18)]);
                    end
                    fcs_vars(m_i) = pa_90(55);
                case 'high'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_90        = [];
                   for dirix = 1:size(dir_vec,1)
                        pa_90(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+48)],[y_coor,poly_yy(dirix:dirix+48)]);
                   end;
                    fcs_vars(m_i) = pa_90(144);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('mean 1-sided distance to canopy in 90° sector to the North-West', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 9)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = mean(poly_d(55:72));
                case 'high'
                    fcs_vars(m_i) = mean(poly_d(144:192));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('max 1-sided distance to canopy in 90° sector to the North-West', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 22)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = max(poly_d(55:72));
                case 'high'
                    fcs_vars(m_i) = max(poly_d(144:192));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('min 1-sided distance to canopy in 90° sector to the North-West', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 35)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = min(poly_d(55:72));
                case 'high'
                    fcs_vars(m_i) = min(poly_d(144:192));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
            %%% parameters calculated over fix 180° sectors
        elseif strcmpi('total open area of 180° sector to the North', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 52)
            switch angular_res
                case 'low'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_180       = [];
                    for dirix = 1:size(dir_vec,1)
                        pa_180(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+36)],[y_coor,poly_yy(dirix:dirix+36)]);
                    end;
                    fcs_vars(m_i) = pa_180(55);
                case 'high'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                   pa_180       = [];
                   for dirix = 1:size(dir_vec,1)
                       pa_180(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+96)],[y_coor,poly_yy(dirix:dirix+96)]);
                   end;
                    fcs_vars(m_i) = pa_180(144);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('mean 1-sided distance to canopy in 180° sector to the North', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 10)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = mean([poly_d(55:72),poly_d(1:19)]);
                case 'high'
                    fcs_vars(m_i) = mean([poly_d(144:192),poly_d(1:48)]);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('max 1-sided distance to canopy in 180° sector to the North', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 23)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = max([poly_d(55:72),poly_d(1:19)]);
                case 'high'
                    fcs_vars(m_i) = max([poly_d(144:192),poly_d(1:48)]);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('min 1-sided distance to canopy in 180° sector to the North', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 36)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = min([poly_d(55:72),poly_d(1:19)]);
                case 'high'
                    fcs_vars(m_i) = min([poly_d(144:192),poly_d(1:48)]);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('total open area of 180° sector to the East', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 53)
            switch angular_res
                case 'low'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_180       = [];
                    for dirix = 1:size(dir_vec,1)
                        pa_180(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+36)],[y_coor,poly_yy(dirix:dirix+36)]);
                    end;
                    fcs_vars(m_i) = pa_90(1);
                case 'high'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_180       = [];
                   for dirix = 1:size(dir_vec,1)
                       pa_180(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+96)],[y_coor,poly_yy(dirix:dirix+96)]);
                   end;
                    fcs_vars(m_i) = pa_90(1);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('mean 1-sided distance to canopy in 180° sector to the East', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 11)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = mean(poly_d(1:37));
                case 'high'
                    fcs_vars(m_i) = mean(poly_d(1:96));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('max 1-sided distance to canopy in 180° sector to the East', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 24)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = max(poly_d(1:37));
                case 'high'
                    fcs_vars(m_i) = max(poly_d(1:96));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('min 1-sided distance to canopy in 180° sector to the East', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 37)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = min(poly_d(1:37));
                case 'high'
                    fcs_vars(m_i) = min(poly_d(1:96));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('total open area of 180° sector to the South', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 54)
            switch angular_res
                case 'low'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_180       = [];
                    for dirix = 1:size(dir_vec,1)
                        pa_180(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+36)],[y_coor,poly_yy(dirix:dirix+36)]);
                    end;
                    fcs_vars(m_i) = pa_90(19);
                case 'high'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_180       = [];
                   for dirix = 1:size(dir_vec,1)
                       pa_180(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+96)],[y_coor,poly_yy(dirix:dirix+96)]);
                   end;
                    fcs_vars(m_i) = pa_90(48);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('mean 1-sided distance to canopy in 180° sector to the South', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 12)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = mean(poly_d(19:55));
                case 'high'
                    fcs_vars(m_i) = mean(poly_d(48:144));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('max 1-sided distance to canopy in 180° sector to the South', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 25)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = max(poly_d(19:55));
                case 'high'
                    fcs_vars(m_i) = max(poly_d(48:144));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
            
        elseif strcmpi('min 1-sided distance to canopy in 180° sector to the South', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 38)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = min(poly_d(19:55));
                case 'high'
                    fcs_vars(m_i) = min(poly_d(48:144));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('total open area of 180° sector to the West', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 55)
            switch angular_res
                case 'low'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_180       = [];
                    for dirix = 1:size(dir_vec,1)
                        pa_180(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+36)],[y_coor,poly_yy(dirix:dirix+36)]);
                    end;
                    fcs_vars(m_i) = pa_90(37);
                case 'high'
                    poly_xx      = [poly_x poly_x];
                    poly_yy      = [poly_y poly_y];
                    pa_180       = [];
                   for dirix = 1:size(dir_vec,1)
                       pa_180(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+96)],[y_coor,poly_yy(dirix:dirix+96)]);
                   end;
                    fcs_vars(m_i) = pa_90(96);
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('mean 1-sided distance to canopy in 180° sector to the West', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 13)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = mean(poly_d(37:72));
                case 'high'
                    fcs_vars(m_i) = mean(poly_d(96:192));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('max 1-sided distance to canopy in 180° sector to the West', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}== 26)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = max(poly_d(37:72));
                case 'high'
                    fcs_vars(m_i) = max(poly_d(96:192));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
        elseif strcmpi('min 1-sided distance to canopy in 180° sector to the West', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}==39)
            switch angular_res
                case 'low'
                    fcs_vars(m_i) = min(poly_d(37:72));
                case 'high'
                    fcs_vars(m_i) = min(poly_d(96:192));
                otherwise
                    error('Error, unknown choice as to angular resolution')
            end
            
            %% > calculate forest canopy structure variables / part 2 (tree height)
            
            %%% identify evaluation domain for additional metric "mean | max tree height"
            
        elseif strcmpi('mean tree height within local evaluation domain', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}==59)
            switch circ_domain
                case 0
                    evix         = find(thr.xm > x_coor - ev_tree_height & thr.xm < x_coor + ev_tree_height & thr.ym > y_coor - ev_tree_height & thr.ym < y_coor + ev_tree_height);
                    fcs_vars(m_i) = mean(thr.zm(evix));
                case 1
                    evix         = find((thr.xm-x_coor).^2 + (thr.ym-y_coor).^2 < ev_tree_height^2);
                    fcs_vars(m_i) = mean(thr.zm(evix));
                otherwise
                    error('Error, unknown choice as to circular domain')
            end
            
        elseif strcmpi('max tree height within local evaluation domain', METRICS{m_i})==1 || ( isnumeric(METRICS{m_i})== 1  && METRICS{m_i}==60)
            switch circ_domain
                case 0
                    evix         = find(thr.xm > x_coor - ev_tree_height & thr.xm < x_coor + ev_tree_height & thr.ym > y_coor - ev_tree_height & thr.ym < y_coor + ev_tree_height);
                    fcs_vars(m_i) = max(thr.zm(evix));
                case 1
                    evix         = find((thr.xm-x_coor).^2 + (thr.ym-y_coor).^2 < ev_tree_height^2);
                    fcs_vars(m_i) = max(thr.zm(evix));
                otherwise
                    error('Error, unknown choice as to circular domain')
            end
        end
        
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% > start arranging forest canopy structure variables
    
    switch meta_raw_data_save
        case 0
            L2CM.vars(:,coordidx)       = fcs_vars;                                % parameter values
            
        case 1
            L2CM.vars(:,coordidx)       = fcs_vars;                                % parameter values
            L2CM.xcoor(coordidx)        = x_coor;                                  % easting of evaluation points
            L2CM.ycoor(coordidx)        = y_coor;                                  % northing of evaluation points
            
        otherwise
            error('Error, unknown choice of "meta_raw_data_save" ')
    end
    
end
clear m_i

%% > finish arranging forest canopy structure variables

switch meta_raw_data_save
    case 0
        L2CM.desc                   = METRICS';
    case 1
        L2CM.desc                   = METRICS';                                % parameter description
        L2CM.file.path              = basefolder;                              % files used
        L2CM.file.LiDAR             = DSM;
        L2CM.file.DTM               = dtm_name;
        L2CM.set.eval_max_dist      = ev_offs;                                 % settings used
        L2CM.set.eval_min_dist      = negl_dist;
        L2CM.set.eval_loc_dist      = ev_tree_height;
        if circ_domain
            L2CM.set.eval_shape       = 'circular';
        else
            L2CM.set.eval_shape       = 'rectangular';
        end
        L2CM.set.height_cutoff      = hg_cutoff;
        L2CM.set.raster_cellsize    = gr_size;
        L2CM.set.raster_smoothing   = smooth_thr;
        L2CM.set.angular_resolution = angular_res;
        L2CM.version                = sversion;
        L2CM.time_minutes           = num2str(toc/60);
        
        % > save forest canopy structure variables
        save([basefolder '\' output_folder '\L2CM_' datestr(now,'yymmddTHHMM') '.mat'],'L2CM','-mat');
        
    otherwise
        error('Error, unknown choice of "meta_raw_data_save" ')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get vector of interpolated values for entire initial dtm and then convert into a matrix (interpolate from evaluation points to the enitre initial domain and initial cellsize)
%  then save as .txtfiles in ARC ASCII format where the file name is based
%  upon the metric name, northing, easting, size of initial domain, and
%  initial cell size

precision  = '%0.6f';
for m_i = 1:length(METRICS)
    clear si interpolated_values
    if sample_interval ~= cellsize
        if  strcmpi (interpolation_method, '')==1  %% if there is no interpolation method set then use 'linear' as default
            interpolation_method = 'linear';
        end
        si                         = scatteredInterpolant(ev_coors(:,1), ev_coors(:,2), (L2CM.vars(m_i,:))',interpolation_method);
        interpolated_values        = si(dtm.x, dtm.y);                 % vector of interpolated CC values
        matrix                     = flipud(reshape(interpolated_values, [ncols,nrows])'); % matrix of interpolated CC values
        
    elseif sample_interval == cellsize
        
        matrix                     = flipud(reshape((L2CM.vars(m_i,:))', [ncols,nrows])');
        
    end
    
    
    if isnumeric(METRICS{m_i})== 1
        text_name =number2name(METRICS{m_i});
        name  = [basefolder '\' output_folder  '/metric_'  num2str(text_name) '_llx_' num2str(easting(1)), '_lly_' num2str(northing(1)) '_nrow_' num2str(nrows) '_ncol_' num2str(ncols), '_cellsize_' num2str(cellsize) '.asc']; % name file
    elseif isnumeric(METRICS{m_i})== 0
        name  = [basefolder  '\' output_folder '/' METRICS{m_i} '_llx_' num2str(easting(1)), '_lly_' num2str(northing(1)) '_nrow_' num2str(nrows) '_ncol_' num2str(ncols), '_cellsize_' num2str(cellsize) '.asc']; % name file
    end
    
    
    fid      =fopen(name,'wt');
    %write header
    fprintf(fid,'%s\t','NCOLS');
    fprintf(fid,'%d\n',ncols);
    
    fprintf(fid,'%s\t','NROWS');
    fprintf(fid,'%d\n',nrows);
    
    fprintf(fid,'%s\t','XLLCORNER');
    fprintf(fid,[precision,'\n'],easting(1));
    
    fprintf(fid,'%s\t','YLLCORNER');
    fprintf(fid,[precision,'\n'],northing(1));
    
    fprintf(fid,'%s\t','CELLSIZE');
    fprintf(fid,[precision,'\n'],cellsize);
    
    fprintf(fid,'%s\t','NODATA_VALUE');
    fprintf(fid,[precision,'\n'],nodata);
    
    fclose(fid);
    dlmwrite(name, matrix,'delimiter',' ','-append');
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%      END ALL      %%%%%%      END ALL      %%%%%%      END ALL      %%%%%%      END ALL      %%%%%%      END ALL      %%%%%%      END ALL      %%%%%%      END ALL      %%%%%%      END ALL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

