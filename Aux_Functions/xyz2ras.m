function [ras] = xyz2ras(fname); 
%  function [ras] = xyz2ras(fname);
%  converts xyz (ASCII) grid to ras grid format
 
% Felix Morsdorf, RSL Zurich, 2006  
  if ~isstruct(fname) 
      [hdr,dat] = hdrload(fname);
      %dat = load(fname);dat(:,1) = [];
  else
    dat = [];
    for i = 1:length(fname)
      disp(['loading x,y,z data from ',fname(i).name])
      [hdr,ndat] = hdrload(fname(i).name);
      ndat = single(ndat);
      dat = single([dat;ndat]);
    end
  end

  % assign variables
  x = single(dat(:,1));y = single(dat(:,2));z = single(dat(:,3));
  clear dat  
  
% determine x-spacing and dimensions

  diffx = diff(x);
  diffx(diffx==0)=[];
  dx = abs(median(diffx));
  xg = min(x):dx:max(x);
  dimx = length(xg);
  %dimx = (max(abs(diffx))/dx)+1;
  
% determine y-spacing and dimensions
  
  diffy = diff(y);
  diffy(diffy==0)=[];  
  dy = abs(median(diffy));
  yg = min(y):dy:max(y);
  dimy = length(yg);
  %dimy = length(z)/dimx;
  % set up new vectors and reshape grid  
  if 0%(fix(dimy)-dimy) == 0 & ~exist('ndat') % no missing values and no mosaic
    ras.z = reshape(z,dimx,dimy)';
    ras.x = xg;
    ras.y = yg;
  else % missing values or mosaicing mode use sort (slow!)
    % TJ: removed warning below
    % disp('Warning : Missing values in Grid, using sort algorithm !');
    rx = xg;
    ry = yg;
    [NX,NY] = meshgrid(rx,ry);
    [m,n] = size(NX);
    NZ = ones(size(NX))*NaN;
    len = length(x);
    for i = 1:len
        ii = x(i) == rx;
        jj = y(i) == ry;
        if sum(ii) == 1 & sum(jj) == 1
          NZ(jj,ii) = z(i);
        end
    end
    ras.z = NZ;
    ras.x = rx;
    ras.y = ry;
  end
  