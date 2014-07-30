% **************************************************************************************
% FMRI_ICA_CLASSIFY: fmri_features4.m: computation of features based on the IC maps and IC timecourses.
% 
% Copyright (C) 2005 - 2007 Jussi Tohka, Institute of Signal Processing, Tampere University of
% Technology, P.O. Box 553, FIN-33101, Finland
% E-mail: jussi.tohka@tut.fi
%
% The method is described in:
% J. Tohka , K. Foerde, A.R. Aron, S. M. Tom, A.W. Toga, and
% R.A. Poldrack. 
% Automatic Independent Component Labeling for Artifact Removal in fMRI.
% NeuroImage , In press , 2007. http://dx.doi.org/10.1016/j.neuroimage.2007.10.013
% ----------------------------------------------------------------------------
% Permission to use, copy, modify, and distribute this software 
% for any purpose and without fee is hereby
% granted, provided that the above copyright notice appear in all
% copies.  The author and Tampere University of Technology make no representations
% about the suitability of this software for any purpose.  It is
% provided "as is" without express or implied warranty.
% ----------------------------------------------------------------------------
% *************************************************************************************
% Computes the features for each IC component (component map, time course + PSD derived from that time course)
% Input :     img    : component map
%             sigft  : component PSD
%             sig    : component timecourse
%             mask   : brain mask
%             bsind  : indeces of voxels in the brain boundary
%             brainind : indeces of voxels well inside the brain
%             varargin = blocked, defaults to 1  
%             flimits 
% Output : feature vectors for each component.
% *************************************************************************************
% PLEASE NOTE that feature numbers in this function are different from
% those in the paper. (The reasons are historical.)
% feature(1) ~ f_1 in the paper
% feature(2) ~ f_4 in the paper
% feature(3) ~ f_3 in the paper
% feature(4) ~ f_2 in the paper
% feature(5) ~ f_5 in the paper
% feature(6) ~ f_6 in the paper
% ************************************************************************************

function feature = fmri_features4(img,sigft,sig,mask,bsind,brainind,varargin);

if(length(varargin) > 0) 
  blocked_design = varargin{1};
else
  blocked_design = 1;
end

if length(varargin) > 1
  flimits = varargin{2};
else
  flimits = [3 5 45]
end

if(blocked_design == 1)
  targetfreq = sum(sigft(9:11));
  lowfreq = sum(sigft(1:flimits(1)));
  reflen = 90;
elseif(blocked_design == 0)
  targetfreq = sum(sigft(flimits(2):flimits(3)));
  lowfreq = sum(sigft(1:flimits(1)));
  reflen = 119;
else 
  targetfreq = sum(sigft((blocked_design - 1):(blocked_design + 1)));
  lowfreq = sum(sigft(1:flimits(1)));
  reflen = 90;
end

sumfreq = sum(sigft);
sz = size(sigft);
% FT
feature(1) = targetfreq/(targetfreq   + lowfreq); 
if(blocked_design > 0)
  feature(1) = ((sz(2)/reflen)^(1/3))*feature(1);
end
% the multiplier tries to compensate the differences of the window length
% (in Hetrzs) for the data that has a different number of time points than
% the training data. It is a heuristic, and it is not needed for the
% event-related case because the window always has approximately the same
% length (in Hertzs). For blocked design data, the window length is 3
% cycles, thus the length of this window depends on the number of time
% points in time series.

sz = size(img);
for i = 1:sz(3)
  pl_voxels(i) = sum(sum(mask(:,:,i)));
end
max_pl_voxels = max(pl_voxels);
pl_ind = find(pl_voxels > 0.6*max_pl_voxels);
pl_ind = min(pl_ind):(max(pl_ind) - mod(max(pl_ind) - min(pl_ind) - 1,2));
  

for i = pl_ind
  maskpl = mask(:,:,i);
  maskplind = find(maskpl(:) > 0);
  imgpl = img(:,:,i);
  imgpl = imgpl(:);  
  planeindex(i) = var(imgpl(maskplind));
end 
s1 = sum(planeindex(pl_ind(1):2:pl_ind(length(pl_ind) - 1)));
s2 = sum(planeindex(pl_ind(2):2:pl_ind(length(pl_ind))));

feature(2) = abs(s1 - s2)/(s1 + s2);
  
maskind = find(mask(:) > 0);

img = img(:);
m = median(img(maskind));
v1 = var(img(maskind));
v2 = var(img(bsind));

feature(3) = (v1 - v2)/v1;

feature(4) = targetfreq/sumfreq;
if(blocked_design > 0)
  feature(4) = ((sz(2)/reflen)^(1/3))*feature(4);
end

[tmp, maxind] = max(abs(diff(sig)));
diffsig = diff(sig);

feature(5) = max(abs(diff(sig)))/mean(abs(diffsig([(1:(maxind - 2)) (maxind + 2):length(diffsig)])));


m = mean(sig);
sig = sig - m;
feature(6) = mean(sig(1:(length(sig) - 1)).*sig(2:length(sig)))/var(sig);

feature(2) = -feature(2);
feature(5) = -feature(5);



