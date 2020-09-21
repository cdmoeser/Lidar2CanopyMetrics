function Lidar2Hemi

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% DESCRIPTION
  %   This scripts creates synthethic hemispherical images from high
  %   resolution LiDAR point cloud data of forest canopy. The script
  %   relates to research published as
  %
  %     D. Moeser, J. Roubinek, P. Schleppi, F. Morsdorf, T. Jonas;
  %     Canopy closure, LAI and radiation transfer from airborne LiDAR
  %     synthetic images; 2014; Agricultural and Forest Meteorology;
  %     doi: 10.1016/j.agrformet.2014.06.008
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
  %   additional toolbox is needed. Add folder Lidar2Hemi including all
  %   subfolders to the Matlab path (to allow access to auxiliarys 
  %   functions called from this script).
  % 
  % IMPLEMENTATION
  %   by David Moeser and Tobias Jonas
  %   WSL Institute for Snow and Avalanche Research SLF
  %   Davos, Switzerland 
  %
  % VERSION / LAST CHANGES
  %   v1.0 / 2015-02-05
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
  %     Evaluation points are coordinates for which synthethic
  %     hemispherical images are to be calculated. These points are defined
  %     from a list of point data. This list must consist of 2 columns
  %     for easting and northing and be saved in a text format.
  %   Note: all three above datasets must feature data in the same
  %     coordinate system (true coordinates). The demonstration dataset
  %     features data in the SwissGrid coordinate syste CH1903/LV03.
  %
  % OUTPUT
  %   this script will generate a synthethic hemispherical image per
  %   coordinate which can be specified in a text file (see user settings).
  %   The images are saved as .gif files in a selectable output file
  %
  % USAGE
  %   adapt user settings as necessary
  %   run Lidar2Hemi from the command line of Matlab
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% USER SETTINGS
  basefolder      = 'C:\scripts\Lidar2Hemi';
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
  
  points          = 'EvPoint_Definitions\laret_low_5_point.txt';
  % relative path\file to text file that contains coordinates for which
  % synthethic hemispherical images are to be calculated. The absolute path
  % to the text file is derived as fullfile(basefolder,points)
  
  output_file     = 'Output_Images';
  % relative path to folder into  which output results will be saved. The
  % absolute path to the output folder is derived as 
  % fullfile(basefolder,output_file)
  
  ev_offs         = 100;
  % this parameter defines the evaluation domain around a point for which 
  % synthethic hemispherical images are to be calculated. Only LiDAR points
  % within this evaluation domain will be considered when calculating a
  % synthethic hemispherical image. The evaluation domain is defined as a
  % diameter around a point, the diameter is to be given in the same unit
  % as the LiDAR data (m in case of the demonstration dataset). Small
  % evaluation domains will allow for quick calculation, but decrease the
  % accuracy of the hemispherical images at high zenith angles (close to
  % to the horizon). Large evaluation domains on the other hand will allow
  % accurate represenation of far distance canopy elements in the 
  % hemispherical images but require more CPU time.
  
  hg_cutoff       = 1.25;
  % this parameter defines a height of canopy elements above the forest
  % floor underneath which all LiDAR points are neglected. A height cutoff
  % is necessary to a) simulate a hemispherical image taken from a certain
  % height above the forst floor, and b) to avoid unrealiable output images
  % due to insufficient point cloud densities as you approach the forest
  % floor (assuming LiDAR data to be taken from above the canopy). 
  
  marker_size     = [5 0.1];
  % this parameter is an import metric that defines the plotting size of  
  % LiDAR points as a function of distance to the camera position. In the 
  % synthetic hemispherical image, a far canopy element should be plotted
  % as a small item, whereas a near canopy element should be represented as
  % a large item. The appearance of the resulting image is significantly 
  % influenced by this parameter. The first number is the plotting size of
  % a single LiDAR point at zero distance to the camera, the second number 
  % is the plotting size of a single LiDAR point at maximum distance to the
  % camera (i.e. ev_offs). Between these two values the plotting size is a
  % linear function of the distance to the camera. Note: an optimal value
  % for this parameter may depend on the point density of the LiDAR data
  % used as well as the output image size. We suggest to adapt this 
  % parameter to allow calculation of realistic LAI values when inputing
  % the resulting synthethic hemispherical images into respective software
  % packages like Hemisfer (cf. figure 4 of Moeser et al., 2014)
  
  radius          = 300;
  % this parameter relates to the output size of the resulting synthethic
  % hemispherical images. The value defines the radius of the 180° image 
  % circle in pixels. Make sure your display has a height of at least 2.5 
  % times the radius.
  
  hidegrid        = 1;                                         
  % default value is 1. Set to 0 in order to display a polar coordinate
  % grid on top of the resulting hemispherical image.
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% CALCULATIONS
  %% > preparatory steps
  
  %%% gather information about computer screen used
  scrsz               = get(0,'ScreenSize');
  screenratio         = scrsz(3) / scrsz(4);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% > read data
  
  %%% read coordinate of points for which synthethic hemispherical images are to be calculated
  coordinates         = dlmread(fullfile(basefolder,points),'');
  coordinate_database = [coordinates(:,1) coordinates(:,2)];               % northing is column 1 and easting is column 2
 
  %%% read LiDAR data in las format
  [raw,hdr]           = readlas(fullfile(basefolder,LiDAR));               % function readlas.m: coutesy of Dr. Felix Morsdorf, raw is the data output, hdr - realtes to the header data
  ras                 = xyz2ras(fullfile(basefolder,DTM));                 % reads DTM data (xyz point cloud) and converts them to ras grid format
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% > preprocess LiDAR and DTM data

  %%% fill gaps in DTM data (ras) 
  [ras.xm,ras.ym]     = meshgrid(ras.x,ras.y);
  ras.xm              = double(ras.xm);                                    % ensure correct format
  ras.ym              = double(ras.ym);                                    % ensure correct format
  ras.xv              = double(reshape(ras.xm,numel(ras.xm),1));           % create vector and ensure correct format
  ras.yv              = double(reshape(ras.ym,numel(ras.ym),1));           % create vector and ensure correct format
  ras.zv              = double(reshape(ras.z ,numel(ras.z ),1));           % create vector and ensure correct format
  nisix               = find(~isnan(ras.zv));                              % find where there are not holes  - this index is used as an input for the below function for a series of three vectors (x,y,z) without holes
  ras.zm              = griddata(ras.xv(nisix),ras.yv(nisix),ras.zv(nisix),ras.xm,ras.ym); % uses above vector set and gives grid bounds (ras.xm, ras.ym) and then linearly interpolates where there is no data

  %%% convert raw to canopy height => raw.z = raw.z - dtm (if not done already)
  if mode(raw.z) > 100
    raw.z = raw.z - interp2(ras.xm,ras.ym,ras.zm,raw.x,raw.y);
  end

  %%% throw away z values outside dtm and hth (if raw.z is less than height cutoff then throw out value)
  raw.z(raw.z<=hg_cutoff) = 0;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% > create synthetic hemispherical image
  
  for coordidx = 1:length(coordinate_database(:,1));                       % loop through coordinates  
    
    %%% create empty output figure which is the size of your screen
    fh = figure('position', [0 0 scrsz(3)  scrsz(4)]);
    figpos = get(fh,'position');
    
    %%% create empty axis with correct size to results in a hemispherical image of radius as specified 
    ah = axes;
    midpoint = [((((scrsz(4)-figpos(4))/2)*screenratio)+figpos(3))/2 (((scrsz(4)-figpos(4))/2)+figpos(4))/2];
    set(ah,'units','pixel');
    set(ah,'position',[midpoint(1)-radius, midpoint(2)-radius*1.15, 2*radius, 2*radius*1.15]);
        
    %%% transfer the global coordinate system into a project ccoordinate system relative to the point in question   
    x_coord(coordidx) = coordinate_database(coordidx,2);
    y_coord(coordidx) = coordinate_database(coordidx,1);
    x_data(:,1)       = raw.x-x_coord(coordidx);
    y_data(:,1)       = raw.y-y_coord(coordidx);

    %%% convert data to a spherical system
    %%% where thet is the counterclockwise angle in the xy plane measured from the positive x axis
    %%% where phi is the elevation angle from the xy plane
    %%% where r is the radius
    [thet(:,1),phi(:,1),r(:,1)] = cart2sph(x_data(:,1),y_data(:,1),raw.z);
    
    %%% flip theta (convert viewing perspective from top-down to bottom-up) 
    thetplot = -thet(:,1);
    
    %%% convert elevation angle phi into zenith angle
    phizplane(:,1) = (pi/2)-phi(:,1);
      
    %%% convert unit of zenith angle to degree and eliminate all values greater than viewing angle (90 degrees)
    phizplanedeg(:,1) = phizplane(:,1)*(180/pi);
    phizplanedeg(phizplanedeg(:,1)>=90,1) = NaN;
      
    %%% identify relevant data, i.e. data within a horizontal radius as specified (ev_offs)
    fix1 = find(r(:,1).*sin(phizplane(:,1)) < ev_offs);

    %%% plot synthetic hemispherical image in case a constant plotting size of LiDAR points is selected (marker_size)
    if length(unique(marker_size)) == 1 
      
      %%% plot all points at once
      ph = polar_special(ah,thetplot(fix1),phizplanedeg(fix1,1),marker_size(1),[0 0 0],hidegrid);
        
    %%% plot synthetic hemispherical image in case a distance-dependant plotting size of LiDAR points is selected (marker_size)
    else
      
      %%% allow 10 categories of marker sizes 
      no_categories = 10;                                                  % 10 provides a good balance between plotting accuracy / time, however this value can be increased if needed
      sizecatbnd = 0:ev_offs/no_categories:ev_offs;
      
      %%% loop through marker size category
      for sizecatix = 1:length(sizecatbnd)-1
        
        %%% identify data in distance categories according to marker size categories
        fix2 = find(r(:,1) >= sizecatbnd(sizecatix) & r(:,1) < sizecatbnd(sizecatix+1)); % select according to 3-d radius (i.e. spherical polar coordinates)
    
        %%% plot only points identifies above
        ph = polar_special(ah,thetplot(fix2),phizplanedeg(fix2,1),(marker_size(1)+(diff(marker_size))*(sizecatix-1)/(no_categories-1)),[0 0 0],hidegrid);
        hold on;
        
      end;
    end;
      
    %%% capture resulting image by recording a screenshot
    scrshot = getframe(fh,[midpoint(1)-radius, midpoint(2)-radius, 2*radius, 2*radius]);
    imwrite(scrshot.cdata(:,:,1),fullfile(basefolder,output_file,['SHI_'  sprintf('%03d',coordidx) '_.gif']),'gif');
    
    %%% close figure and clear workspace to allow clear variables for next calculation
    close
    clear fix1 fix2 x_coord y_coord x_data y_data thet thetplot phi r phizplane phizplanedeg
      
  end


    