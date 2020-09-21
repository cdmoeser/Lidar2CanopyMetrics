function [raw] = subsetraw(raw,ii);
%  function [raw] = subsetraw(raw,ii);
%  Subset structure raw with logical array ii of length(raw.x)
%  works for all field of raw
  
% Felix Morsdorf, RSL Zurich, April 2009
  
  if islogical(ii) 
    fnames = fieldnames(raw);
    for i = 1:length(fnames)
      eval(['raw.',fnames{i},' = raw.',fnames{i},'(ii);']);
      eval(['raw.',fnames{i},' = raw.',fnames{i},'(:);']);
    end
  else
    fnames = fieldnames(raw);
    for i = 1:length(fnames)
      eval(['raw.',fnames{i},' = raw.',fnames{i},'(ii);']);
    end
  end
  
  