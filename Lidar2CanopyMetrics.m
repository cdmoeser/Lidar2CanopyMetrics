function Lidar2CanopyMetrics

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
  %% USER SETTINGS
  basefolder      = 'C:\scripts\Lidar2CanopyMetrics';
  % path to folder which contains this script as well as all subfolders 
  % required to run this script 
  
  LiDAR           = 'LiDAR_Data\laret_low.las';
  % relative path\file to data file that contains LiDAR data (.las format).
  % The absolute path to the LiDAR data file is derived as 
  % fullfile(basefolder,LiDAR)
  
  DTM             = 'DTM_Data\laret_low.txt';
  % relative path\file to data file that contains DTM data (.txt format).
  % The absolute path to the DTM data file is derived as 
  % fullfile(basefolder,DTM)
  
  points          = 'EvPoint_Definitions\laret_low_16_point.txt';
  % relative path\file to text file that contains coordinates for which
  % forest canopy structure metrics are to be calculated. The absolute path
  % to the text file is derived as fullfile(basefolder,points)
  
  output_file     = 'Output_Data';
  % relative path to folder into  which output results will be saved. The
  % absolute path to the output folder is derived as 
  % fullfile(basefolder,output_file)
  
  ev_offs         = 300;
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
  % the additional metric "average tree height" is to be calculated. Only
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
  % set to 1 if you wish to display a plot per evaluation point, set to 0 
  % otherwise. The plot will show an analysis of the gap space around 
  % each evaluation point.
  
  sversion        = 'v1.2 / 2016-05-30';
  % string, will be save in the results file for reference purposes 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% CALCULATIONS
  %% > preparatory steps

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% > read data
  
  %%% read coordinate of evaluation points for which synthethic hemispherical images are to be calculated
  coordinates           = dlmread(fullfile(basefolder,points),'');
  coordinate_database   = [coordinates(:,1) coordinates(:,2)];             % northing is column 1 and easting is column 2
 
  %%% read LiDAR data in las format
  [raw,hdr]             = readlas(fullfile(basefolder,LiDAR));             % function readlas.m: coutesy of Dr. Felix Morsdorf, raw is the data output, hdr - realtes to the header data
  ras                   = xyz2ras(fullfile(basefolder,DTM));               % reads DTM data (xyz point cloud) and converts them to ras grid format
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% > preprocess LiDAR and DTM data

  %%% fill gaps in DTM data (ras) 
  [ras.xm,ras.ym]       = meshgrid(ras.x,ras.y);
  ras.xm                = double(ras.xm);                                  % ensure correct format
  ras.ym                = double(ras.ym);                                  % ensure correct format
  ras.xv                = double(reshape(ras.xm,numel(ras.xm),1));         % create vector and ensure correct format
  ras.yv                = double(reshape(ras.ym,numel(ras.ym),1));         % create vector and ensure correct format
  ras.zv                = double(reshape(ras.z ,numel(ras.z ),1));         % create vector and ensure correct format
  nisix                 = find(~isnan(ras.zv));                            % find where there are not holes  - this index is used as an input for the below function for a series of three vectors (x,y,z) without holes
  ras.zm                = griddata(ras.xv(nisix),ras.yv(nisix),ras.zv(nisix),ras.xm,ras.ym); % uses above vector set and gives grid bounds (ras.xm, ras.ym) and then linearly interpolates where there is no data

  %%% convert raw to canopy height => raw.z = raw.z - dtm (if not done already)
  if mode(raw.z) > 100
    raw.z = raw.z - interp2(ras.xm,ras.ym,ras.zm,raw.x,raw.y);
  end
  
  %%% remove outliers, dismiss canopy elements with unreasonable heights (< 0m / > 100m)
  ii                    = raw.z > 0 & raw.z < 100;
  raw                   = subsetraw(raw,ii);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% > prepare binarized canopy height raster
  
  %%% create canopy height raster for evaluation domain
  x_coors               = coordinate_database(:,2);
  y_coors               = coordinate_database(:,1);
  xlim                  = [min(x_coors)-ev_offs,max(x_coors)+ev_offs];     % defines rectangular evaluation domain
  ylim                  = [min(y_coors)-ev_offs,max(y_coors)+ev_offs];     % defines rectangular evaluation domain
  [thr.xm,thr.ym]       = meshgrid(xlim(1):gr_size:xlim(2),ylim(1):gr_size:ylim(2)); % defines grid coordinates
  nisix                 = find(~isnan(raw.z));                               
  thr.zm                = griddata(raw.x(nisix),raw.y(nisix),raw.z(nisix),thr.xm,thr.ym);
  
  %%% fill all nan cells in the canopy height raster with 0, except for entire rows/columns with nan values
  nanix                 = isnan(thr.zm);                                   % identify nan cells
  allnancols            = find(all(nanix,1));                              % identify columns with all values being nan
  allnanrows            = find(all(nanix,2));                              % identify rows with all values being nan
  thr.zm(isnan(thr.zm)) = 0;                                               % fill isolated gaps with 0 (no canopy element)
  thr.zm(allnanrows,:)  = NaN;                                             % re-establish nan values for entire rows
  thr.zm(:,allnancols)  = NaN;                                             % re-establish nan values for entire columns
    
  %%% binarize canopy height raster
  thr.zb                = thr.zm;                                            
  thr.zb(thr.zm < hg_cutoff) = 0;
  thr.zb(thr.zm >= hg_cutoff) = 1;
    
  %%% smooth binarized canopy height raster using a moving 3 × 3 neighborhood cell mean filter
  thr.zb = smooth2a(thr.zb,1,1);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% > prepare calculating forest canopy structure metrics
  
  for coordidx = 1:length(coordinate_database(:,1));                       % loop through coordinates of evaluation points
    
    %%% create canopy height raster for evaluation domain
    x_coor            = coordinate_database(coordidx,2);
    y_coor            = coordinate_database(coordidx,1);
    
    %%% identify position of evaluation point within canopy height raster 
    nox               = length(thr.xm);                                    
    noy               = length(thr.ym);                                     
    xll               = max(thr.xm(1,thr.xm(1,1:end) <= x_coor)); 
    yll               = max(thr.ym(thr.ym(1:end,1) <= y_coor,1));
    if isempty(xll) || x_coor-xll > gr_size || isempty(yll) || y_coor-yll > gr_size
      error('Error locating position of point for which forest canopy structure metrics are to be calculated within canopy height raster, this should not happen')
    else
      xix = find(abs(thr.xm(1,:)-xll)<eps);
      yix = find(abs(thr.ym(:,1)-yll)<eps);
    end;
    negl_ix = round(negl_dist/gr_size);
    
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
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% > analyze gap space around evaluation point
    
    poly_x = [];
    poly_y = [];
    poly_d = [];
    
    for dirix = 1:size(dir_vec,1)                                          % loop through directions [N, NNE, NE, ENE, E, ESE, SE, SSE, S, SSW, SW, WSW, W, WNW, NW, NNW
      
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
          end;
        end;
      end;
    end;
    
    %%% optionally restrict results to circular domain
    if circ_domain
      for dirix = 1:size(dir_vec,1)
        if poly_d(dirix) > ev_offs
          poly_x(dirix) = x_coor+(poly_x(dirix)-x_coor)*ev_offs/poly_d(dirix);
          poly_y(dirix) = y_coor+(poly_y(dirix)-y_coor)*ev_offs/poly_d(dirix);
          poly_d(dirix) = ev_offs;
        end;
      end;
    end;
    
    %%% create a plot with an analysis of the gap space around the evaluation point
    if doplot
      fh = figure;
      ph = pcolor(thr.xm,thr.ym,thr.zb); hold all
      set(ph,'linestyle','none');
      plot(x_coor, y_coor,'.K');
      pg = fill(poly_x+1/2*gr_size,poly_y+1/2*gr_size,'w');
      set(pg,'EdgeColor',[0 1 1],'FaceColor','none','LineWidth',2)
      for dirix = 1:size(dir_vec,1)
        plot([x_coor,poly_x(dirix)]+1/2*gr_size,[y_coor,poly_y(dirix)]+1/2*gr_size,'g-');
      end;
      set(gca,'xlim',[min(thr.xm(~isnan(thr.zb))),max(thr.xm(~isnan(thr.zb)))],'ylim',[min(thr.ym(~isnan(thr.zb))),max(thr.ym(~isnan(thr.zb)))]);
      drawnow;
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% > calculate forest canopy structure variables / part 1 (vector search)
    
    fcs_vars = [];
    fcs_desc = {};
    switch angular_res
    case 'low' % search for canopy elements every 5° (72 directions)
      
      %%% parameters calculated over entire visisble gap space
      fcs_desc{1}  = 'average 2-sided distance to canopy';
      fcs_vars(1)  = mean(poly_d(1:36)+poly_d(37:72));
      fcs_desc{2}  = 'maximum 2-sided distance to canopy';
      fcs_vars(2)  = max(poly_d(1:36)+poly_d(37:72));
      fcs_desc{3}  = 'minimum 2-sided distance to canopy';
      fcs_vars(3)  = min(poly_d(1:36)+poly_d(37:72));
      fcs_desc{4}  = 'average 1-sided distance to canopy';
      fcs_vars(4)  = mean(poly_d(1:72));
      fcs_desc{5}  = 'maximum 1-sided distance to canopy';
      fcs_vars(5)  = max(poly_d(1:72));
      fcs_desc{6}  = 'minimum 1-sided distance to canopy';
      fcs_vars(6)  = min(poly_d(1:72));
      fcs_desc{7}  = 'total visible open area';
      fcs_vars(7)  = polyarea(poly_x,poly_y);
      
      %%% parameters calculated over moving angular sectors
      poly_xx      = [poly_x poly_x];
      poly_yy      = [poly_y poly_y];
      pa_45        = [];
      for dirix = 1:size(dir_vec,1)
        pa_45(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+9)],[y_coor,poly_yy(dirix:dirix+9)]);
      end;
      fcs_desc{8}  = 'maximum visible open area of any 45 degree sector';
      fcs_vars(8)  = max(pa_45);
      pa_90        = [];
      for dirix = 1:size(dir_vec,1)
        pa_90(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+18)],[y_coor,poly_yy(dirix:dirix+18)]);
      end;
      fcs_desc{9}  = 'maximum visible open area of any 90 degree sector';
      fcs_vars(9)  = max(pa_90);
      pa_180       = [];
      for dirix = 1:size(dir_vec,1)
        pa_180(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+36)],[y_coor,poly_yy(dirix:dirix+36)]);
      end;
      fcs_desc{10}  = 'maximum visible open area of any 180 degree sector';
      fcs_vars(10)  = max(pa_180);
        
      %%% parameters calculated over fix 90° sectors
      fcs_desc{11}  = 'visible open area of 90° sector to the North';
      fcs_vars(11)  = pa_90(64);
      fcs_desc{12} = 'average 1-sided distance to canopy in 90° sector to the North';
      fcs_vars(12) = mean([poly_d(64:72),poly_d(1:10)]);
      fcs_desc{13} = 'max 1-sided distance to canopy in 90° sector to the North';
      fcs_vars(13) = max([poly_d(64:72),poly_d(1:10)]);
      fcs_desc{14} = 'min 1-sided distance to canopy in 90° sector to the North';
      fcs_vars(14) = min([poly_d(64:72),poly_d(1:10)]);
      fcs_desc{15} = 'visible open area of 90° sector to the North-East';
      fcs_vars(15) = pa_90(1);
      fcs_desc{16} = 'average 1-sided distance to canopy in 90° sector to the North-East';
      fcs_vars(16) = mean(poly_d(1:19));
      fcs_desc{17} = 'max 1-sided distance to canopy in 90° sector to the North-East';
      fcs_vars(17) = max(poly_d(1:19));
      fcs_desc{18} = 'min 1-sided distance to canopy in 90° sector to the North-East';
      fcs_vars(18) = min(poly_d(1:19));
      fcs_desc{19} = 'visible open area of 90° sector to the East';
      fcs_vars(19) = pa_90(10);
      fcs_desc{20} = 'average 1-sided distance to canopy in 90° sector to the East';
      fcs_vars(20) = mean(poly_d(10:28));
      fcs_desc{21} = 'max 1-sided distance to canopy in 90° sector to the East';
      fcs_vars(21) = max(poly_d(10:28));
      fcs_desc{22} = 'min 1-sided distance to canopy in 90° sector to the East';
      fcs_vars(22) = min(poly_d(10:28));
      fcs_desc{23} = 'visible open area of 90° sector to the South-East';
      fcs_vars(23) = pa_90(19);
      fcs_desc{24} = 'average 1-sided distance to canopy in 90° sector to the South-East';
      fcs_vars(24) = mean(poly_d(19:37));
      fcs_desc{25} = 'max 1-sided distance to canopy in 90° sector to the South-East';
      fcs_vars(25) = max(poly_d(19:37));
      fcs_desc{26} = 'min 1-sided distance to canopy in 90° sector to the South-East';
      fcs_vars(26) = min(poly_d(19:37));
      fcs_desc{27} = 'visible open area of 90° sector to the South';
      fcs_vars(27) = pa_90(28);
      fcs_desc{28} = 'average 1-sided distance to canopy in 90° sector to the South';
      fcs_vars(28) = mean(poly_d(28:46));
      fcs_desc{29} = 'max 1-sided distance to canopy in 90° sector to the South';
      fcs_vars(29) = max(poly_d(28:46));
      fcs_desc{30} = 'min 1-sided distance to canopy in 90° sector to the South';
      fcs_vars(30) = min(poly_d(28:46));
      fcs_desc{31} = 'visible open area of 90° sector to the South-West';
      fcs_vars(31) = pa_90(37);
      fcs_desc{32} = 'average 1-sided distance to canopy in 90° sector to the South-West';
      fcs_vars(32) = mean(poly_d(37:55));
      fcs_desc{33} = 'max 1-sided distance to canopy in 90° sector to the South-West';
      fcs_vars(33) = max(poly_d(37:55));
      fcs_desc{34} = 'min 1-sided distance to canopy in 90° sector to the South-West';
      fcs_vars(34) = min(poly_d(37:55));
      fcs_desc{35} = 'visible open area of 90° sector to the West';
      fcs_vars(35) = pa_90(46);
      fcs_desc{36} = 'average 1-sided distance to canopy in 90° sector to the West';
      fcs_vars(36) = mean(poly_d(46:64));
      fcs_desc{37} = 'max 1-sided distance to canopy in 90° sector to the West';
      fcs_vars(37) = max(poly_d(46:64));
      fcs_desc{38} = 'min 1-sided distance to canopy in 90° sector to the West';
      fcs_vars(38) = min(poly_d(46:64));
      fcs_desc{39} = 'visible open area of 90° sector to the North-West';
      fcs_vars(39) = pa_90(55);
      fcs_desc{40} = 'average 1-sided distance to canopy in 90° sector to the North-West';
      fcs_vars(40) = mean(poly_d(55:72));
      fcs_desc{41} = 'max 1-sided distance to canopy in 90° sector to the North-West';
      fcs_vars(41) = max(poly_d(55:72));
      fcs_desc{42} = 'min 1-sided distance to canopy in 90° sector to the North-West';
      fcs_vars(42) = min(poly_d(55:72));
      
      %%% parameters calculated over fix 180° sectors
      fcs_desc{43} = 'visible open area of 180° sector to the North';
      fcs_vars(43) = pa_180(55);
      fcs_desc{44} = 'average 1-sided distance to canopy in 180° sector to the North';
      fcs_vars(44) = mean([poly_d(55:72),poly_d(1:19)]);
      fcs_desc{45} = 'max 1-sided distance to canopy in 180° sector to the North';
      fcs_vars(45) = max([poly_d(55:72),poly_d(1:19)]);
      fcs_desc{46} = 'min 1-sided distance to canopy in 180° sector to the North';
      fcs_vars(46) = min([poly_d(55:72),poly_d(1:19)]);
      fcs_desc{47} = 'visible open area of 180° sector to the East';
      fcs_vars(47) = pa_90(1);
      fcs_desc{48} = 'average 1-sided distance to canopy in 180° sector to the East';
      fcs_vars(48) = mean(poly_d(1:37));
      fcs_desc{49} = 'max 1-sided distance to canopy in 180° sector to the East';
      fcs_vars(49) = max(poly_d(1:37));
      fcs_desc{50} = 'min 1-sided distance to canopy in 180° sector to the East';
      fcs_vars(50) = min(poly_d(1:37));
      fcs_desc{51} = 'visible open area of 180° sector to the South';
      fcs_vars(51) = pa_90(19);
      fcs_desc{52} = 'average 1-sided distance to canopy in 180° sector to the South';
      fcs_vars(52) = mean(poly_d(19:55));
      fcs_desc{53} = 'max 1-sided distance to canopy in 180° sector to the South';
      fcs_vars(53) = max(poly_d(19:55));
      fcs_desc{54} = 'min 1-sided distance to canopy in 180° sector to the South';
      fcs_vars(54) = min(poly_d(19:55));
      fcs_desc{55} = 'visible open area of 180° sector to the West';
      fcs_vars(55) = pa_90(37);
      fcs_desc{56} = 'average 1-sided distance to canopy in 180° sector to the West';
      fcs_vars(56) = mean(poly_d(37:72));
      fcs_desc{57} = 'max 1-sided distance to canopy in 180° sector to the West';
      fcs_vars(57) = max(poly_d(37:72));
      fcs_desc{58} = 'min 1-sided distance to canopy in 180° sector to the West';
      fcs_vars(58) = min(poly_d(37:72));
      
    case 'high' % search for canopy elements every 2° (192 directions)
      
      %%% parameters calculated over entire visisble gap space
      fcs_desc{1}  = 'average 2-sided distance to canopy';
      fcs_vars(1)  = mean(poly_d(1:96)+poly_d(97:192));
      fcs_desc{2}  = 'maximum 2-sided distance to canopy';
      fcs_vars(2)  = max(poly_d(1:96)+poly_d(97:192));
      fcs_desc{3}  = 'minimum 2-sided distance to canopy';
      fcs_vars(3)  = min(poly_d(1:96)+poly_d(97:192));
      fcs_desc{4}  = 'average 1-sided distance to canopy';
      fcs_vars(4)  = mean(poly_d(1:192));
      fcs_desc{5}  = 'maximum 1-sided distance to canopy';
      fcs_vars(5)  = max(poly_d(1:192));
      fcs_desc{6}  = 'minimum 1-sided distance to canopy';
      fcs_vars(6)  = min(poly_d(1:192));
      fcs_desc{7}  = 'total visible open area';
      fcs_vars(7)  = polyarea(poly_x,poly_y);
      
      %%% parameters calculated over moving angular sectors
      poly_xx      = [poly_x poly_x];
      poly_yy      = [poly_y poly_y];
      pa_45        = [];
      for dirix = 1:size(dir_vec,1)
        pa_45(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+24)],[y_coor,poly_yy(dirix:dirix+24)]);
      end;
      fcs_desc{8}  = 'maximum visible open area of any 45 degree sector';
      fcs_vars(8)  = max(pa_45);
      pa_90        = [];
      for dirix = 1:size(dir_vec,1)
        pa_90(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+48)],[y_coor,poly_yy(dirix:dirix+48)]);
      end;
      fcs_desc{9}  = 'maximum visible open area of any 90 degree sector';
      fcs_vars(9)  = max(pa_90);
      pa_180       = [];
      for dirix = 1:size(dir_vec,1)
        pa_180(dirix) = polyarea([x_coor,poly_xx(dirix:dirix+96)],[y_coor,poly_yy(dirix:dirix+96)]);
      end;
      fcs_desc{10}  = 'maximum visible open area of any 180 degree sector';
      fcs_vars(10)  = max(pa_180);
      
      %%% parameters calculated over fix 90° sectors
      fcs_desc{11}  = 'visible open area of 90° sector to the North';
      fcs_vars(11)  = pa_90(168);
      fcs_desc{12} = 'average 1-sided distance to canopy in 90° sector to the North';
      fcs_vars(12) = mean([poly_d(168:192),poly_d(1:24)]);
      fcs_desc{13} = 'max 1-sided distance to canopy in 90° sector to the North';
      fcs_vars(13) = max([poly_d(168:192),poly_d(1:24)]);
      fcs_desc{14} = 'min 1-sided distance to canopy in 90° sector to the North';
      fcs_vars(14) = min([poly_d(168:192),poly_d(1:24)]);
      fcs_desc{15} = 'visible open area of 90° sector to the North-East';
      fcs_vars(15) = pa_90(1);
      fcs_desc{16} = 'average 1-sided distance to canopy in 90° sector to the North-East';
      fcs_vars(16) = mean(poly_d(1:48));
      fcs_desc{17} = 'max 1-sided distance to canopy in 90° sector to the North-East';
      fcs_vars(17) = max(poly_d(1:48));
      fcs_desc{18} = 'min 1-sided distance to canopy in 90° sector to the North-East';
      fcs_vars(18) = min(poly_d(1:48));
      fcs_desc{19} = 'visible open area of 90° sector to the East';
      fcs_vars(19) = pa_90(24);
      fcs_desc{20} = 'average 1-sided distance to canopy in 90° sector to the East';
      fcs_vars(20) = mean(poly_d(24:72));
      fcs_desc{21} = 'max 1-sided distance to canopy in 90° sector to the East';
      fcs_vars(21) = max(poly_d(24:72));
      fcs_desc{22} = 'min 1-sided distance to canopy in 90° sector to the East';
      fcs_vars(22) = min(poly_d(24:72));
      fcs_desc{23} = 'visible open area of 90° sector to the South-East';
      fcs_vars(23) = pa_90(48);
      fcs_desc{24} = 'average 1-sided distance to canopy in 90° sector to the South-East';
      fcs_vars(24) = mean(poly_d(48:96));
      fcs_desc{25} = 'max 1-sided distance to canopy in 90° sector to the South-East';
      fcs_vars(25) = max(poly_d(48:96));
      fcs_desc{26} = 'min 1-sided distance to canopy in 90° sector to the South-East';
      fcs_vars(26) = min(poly_d(48:96));
      fcs_desc{27} = 'visible open area of 90° sector to the South';
      fcs_vars(27) = pa_90(72);
      fcs_desc{28} = 'average 1-sided distance to canopy in 90° sector to the South';
      fcs_vars(28) = mean(poly_d(72:120));
      fcs_desc{29} = 'max 1-sided distance to canopy in 90° sector to the South';
      fcs_vars(29) = max(poly_d(72:120));
      fcs_desc{30} = 'min 1-sided distance to canopy in 90° sector to the South';
      fcs_vars(30) = min(poly_d(72:120));
      fcs_desc{31} = 'visible open area of 90° sector to the South-West';
      fcs_vars(31) = pa_90(96);
      fcs_desc{32} = 'average 1-sided distance to canopy in 90° sector to the South-West';
      fcs_vars(32) = mean(poly_d(96:144));
      fcs_desc{33} = 'max 1-sided distance to canopy in 90° sector to the South-West';
      fcs_vars(33) = max(poly_d(96:144));
      fcs_desc{34} = 'min 1-sided distance to canopy in 90° sector to the South-West';
      fcs_vars(34) = min(poly_d(96:144));
      fcs_desc{35} = 'visible open area of 90° sector to the West';
      fcs_vars(35) = pa_90(120);
      fcs_desc{36} = 'average 1-sided distance to canopy in 90° sector to the West';
      fcs_vars(36) = mean(poly_d(120:168));
      fcs_desc{37} = 'max 1-sided distance to canopy in 90° sector to the West';
      fcs_vars(37) = max(poly_d(120:168));
      fcs_desc{38} = 'min 1-sided distance to canopy in 90° sector to the West';
      fcs_vars(38) = min(poly_d(120:168));
      fcs_desc{39} = 'visible open area of 90° sector to the North-West';
      fcs_vars(39) = pa_90(144);
      fcs_desc{40} = 'average 1-sided distance to canopy in 90° sector to the North-West';
      fcs_vars(40) = mean(poly_d(144:192));
      fcs_desc{41} = 'max 1-sided distance to canopy in 90° sector to the North-West';
      fcs_vars(41) = max(poly_d(144:192));
      fcs_desc{42} = 'min 1-sided distance to canopy in 90° sector to the North-West';
      fcs_vars(42) = min(poly_d(144:192));
      
      %%% parameters calculated over fix 180° sectors
      fcs_desc{43} = 'visible open area of 180° sector to the North';
      fcs_vars(43) = pa_180(144);
      fcs_desc{44} = 'average 1-sided distance to canopy in 180° sector to the North';
      fcs_vars(44) = mean([poly_d(144:192),poly_d(1:48)]);
      fcs_desc{45} = 'max 1-sided distance to canopy in 180° sector to the North';
      fcs_vars(45) = max([poly_d(144:192),poly_d(1:48)]);
      fcs_desc{46} = 'min 1-sided distance to canopy in 180° sector to the North';
      fcs_vars(46) = min([poly_d(144:192),poly_d(1:48)]);
      fcs_desc{47} = 'visible open area of 180° sector to the East';
      fcs_vars(47) = pa_90(1);
      fcs_desc{48} = 'average 1-sided distance to canopy in 180° sector to the East';
      fcs_vars(48) = mean(poly_d(1:96));
      fcs_desc{49} = 'max 1-sided distance to canopy in 180° sector to the East';
      fcs_vars(49) = max(poly_d(1:96));
      fcs_desc{50} = 'min 1-sided distance to canopy in 180° sector to the East';
      fcs_vars(50) = min(poly_d(1:96));
      fcs_desc{51} = 'visible open area of 180° sector to the South';
      fcs_vars(51) = pa_90(48);
      fcs_desc{52} = 'average 1-sided distance to canopy in 180° sector to the South';
      fcs_vars(52) = mean(poly_d(48:144));
      fcs_desc{53} = 'max 1-sided distance to canopy in 180° sector to the South';
      fcs_vars(53) = max(poly_d(48:144));
      fcs_desc{54} = 'min 1-sided distance to canopy in 180° sector to the South';
      fcs_vars(54) = min(poly_d(48:144));
      fcs_desc{55} = 'visible open area of 180° sector to the West';
      fcs_vars(55) = pa_90(96);
      fcs_desc{56} = 'average 1-sided distance to canopy in 180° sector to the West';
      fcs_vars(56) = mean(poly_d(96:192));
      fcs_desc{57} = 'max 1-sided distance to canopy in 180° sector to the West';
      fcs_vars(57) = max(poly_d(96:192));
      fcs_desc{58} = 'min 1-sided distance to canopy in 180° sector to the West';
      fcs_vars(58) = min(poly_d(96:192));
      
    otherwise
      error('Error, unknown choice as to angular resolution')
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% > calculate forest canopy structure variables / part 2 (tree height)

    %%% identify evaluation domain for additional metric "average tree height"
    switch circ_domain
    case 0
      evix         = find(thr.xm > x_coor - ev_tree_height & thr.xm < x_coor + ev_tree_height & thr.ym > y_coor - ev_tree_height & thr.ym < y_coor + ev_tree_height);
      fcs_desc{59} = 'average tree height within local evaluation domain';
      fcs_vars(59) = mean(thr.zm(evix));
      fcs_desc{60} = 'maximum tree height within local evaluation domain';
      fcs_vars(60) = max(thr.zm(evix));
    case 1
      evix         = find((thr.xm-x_coor).^2 + (thr.ym-y_coor).^2 < ev_tree_height^2);
      fcs_desc{59} = 'average tree height within local evaluation domain';
      fcs_vars(59) = mean(thr.zm(evix));
      fcs_desc{60} = 'maximum tree height within local evaluation domain';
      fcs_vars(60) = max(thr.zm(evix));
    otherwise
      error('Error, unknown choice as to circular domain')
    end;   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% > arrange forest canopy structure variables
    L2CM.desc                   = fcs_desc';                               % parameter description
    L2CM.vars(:,coordidx)       = fcs_vars;                                % parameter values
    L2CM.xcoor(coordidx)        = x_coor;                                  % easting of evaluation points
    L2CM.ycoor(coordidx)        = y_coor;                                  % northing of evaluation points
    L2CM.file.path              = basefolder;                              % files used
    L2CM.file.LiDAR             = LiDAR;
    L2CM.file.DTM               = DTM;
    L2CM.file.points            = points;
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
        
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% > save forest canopy structure variables
  save(fullfile(basefolder,output_file,['L2CM_' datestr(now,'yymmddTHHMM') '.mat']),'L2CM','-mat');
end  
    

    
    
    
    
    

  
  


 


    