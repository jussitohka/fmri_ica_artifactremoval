% FMRI_ICA_TRAINING Version 1.2.1:
% Trains a global decision tree classifier to assign the ICA components 
% of an fMR image to "signal" and "noise" for artifact removal
% Copyright (C) 2005 - 2009 Jussi Tohka, Institute of Signal Processing, Tampere University of
% Technology, P.O. Box 553, FIN-33101, Finland
% E-mail: jussi.tohka@tut.fi
%
% The method is described in:
% J. Tohka, K. Foerde, A.R. Aron, S. M. Tom, A.W. Toga, and
% R.A. Poldrack. 
% Automatic Independent Component Labeling for Artifact Removal in fMRI.
% NeuroImage, Volume 39, Issue 3, 1 February 2008, Pages 1227-1245, 2008.
%  http://dx.doi.org/10.1016/j.neuroimage.2007.10.013
% ----------------------------------------------------------------------------
% Version 1.1 23rd Jan 2008: The possibility to read gzipped files was added (JT). 
% Version 1.2 27th Oct 2008: 1) The possibility to deal with classes with no 
%                            training data was added (JT).
%                            2) the possibility to define the noise classes in
%                               training data was added
%                            3) minor bug corrections took place
% Version 1.2.1. 19th Nov 2009: Updated the function fmri_readicanii to use 
%               read_untouch_nii instead of read_nii of the previous version. 
%               This reflects changes to the NIFTI-Matlab Toolbox by Jimmy Shen,
%               the current version is up-to-date with the NIFTI Toolbox dated 
  %               20090909. (JT)    
% ----------------------------------------------------------------------------
% Permission to use, copy, modify, and distribute this software 
% for any purpose and without fee is hereby
% granted, provided that the above copyright notice appear in all
% copies.  The author and Tampere University of Technology make no representations
% about the suitability of this software for any purpose.  It is
% provided "as is" without express or implied warranty.
% ----------------------------------------------------------------------------
% Requires functions       : fmri_features4, global_decision_tree_v2
%
% ----------------------------------------------------------------------------
%
% Usage   :  [global_decision_tree_classifier, rejected_clusters] = 
%  
%            fmri_ica_training(training_data_fn,np_thr,np_testlevels,blocked,tr,flimits,input_rejected_clusters);
%          
%                          
% ************************************************************************************
% Input    :    training_data_fn : the name of the text file specifying the training data. 
%                                  See the documentation for the specification of this file.   
%               np_thr           : the desired Neyman Pearson threshold(s). If needed to train classifiers 
%                                  for multiple NP thresholds, give a vector of threshold values
%                                  for which the classifiers then are trained. This does not require extra
%                                  computation time.   
%               np_testlevels    : The NP thresholds for which the element classifiers are trained.
%                                  Defaults to [0.0100 0.0150 0.0200 0.0300 0.0400 0.0500 0.0700 0.0900]. 
%                                  If you are happy with the default, you can give an empty matrix.
%               blocked          : 0 if event related design (target
%                                  frequency assumed to be 0.01 Hz -  0.1 Hz)
%                                : 1 if blocked design with the target frequency of 10 cycles
%                                : n if blocked design with the target frequency of n cycles
%              tr(optional)  : The repetition time (TR) in ms. Defaults to
%                              2000. If you're using a considerably
%                              different TR than 2000, the classifier may
%                              not work well for the event related case. 
%              flimits(optional) : A 3 component vector which sets 1) low
%                                  frequency limit, 2) the lower limit of the target
%                                  frequency and 3) the upper limit of the target
%                                  frequency. All should be given in Hertzs. 2) and 3) are
%                                  ignored for the blocked designs. Give an empty matrix to go with the defaults, 
%                                  if you want to use the rejected_clusters option.  
%    input_rejected_clusters(optional) : If you want to specify the noise classes of the rejected components
%                                  you can give them here. See the output section. 
% ************************************************************************************
%  Output: 
%    global_decision_tree_classifier:  A struct containing information of the training process.
%                                      global_decision_tree_classifier.classifier is a cell array of 
%                                      the trained GDT classifiers with different NP thresholds. 
%                                      You can give one of these as an input to fmri_ica_classify.m
%    rejected_clusters               : gives the noise classes of rejected components for debugging purposes. 
%                                      This is a cell matrix with the number of cells equalling the number 
%                                      of training timeseries. Each cell is a matrix n by 2 matrix whose each 
%                                      row gives the number of the rejected component and its noise class. 
%                                      i.e. rejected_clusters{3}(2,1) gives the number of 2nd rejected ICA component of 
%                                      the 3rd time-series and  rejected_clusters{3}(2,2) gives its noise class as 
%                                      determined by the clustering algorithm described in the Appendix of the paper. 
%                                      This
%                                      can be used to check if you agree with the clustering result. If you 
%                                      don't agree with the clustering you can change the noise classes here, and 
%                                      give it as a final input to the function.       
% ************************************************************************************

function [global_decision_tree_classifier,rejected_clusters] = fmri_ica_training(training_data_fn,np_thr,np_testlevels,blocked,varargin);

% compile the filenames and component classes for the training data

  fid = fopen(training_data_fn,'r');
  n_of_timeseries = 1;
  stmp = fscanf(fid,'%s',1);
  while(~isempty(stmp))
     icadir{n_of_timeseries} = stmp;
     g = fgetl(fid);
     rej{n_of_timeseries} = str2num(g);
     n_of_timeseries = n_of_timeseries + 1;
     stmp = fscanf(fid,'%s',1);
  end
  fclose(fid);
  n_of_timeseries = n_of_timeseries - 1;
% get TR and flimits
  if length(varargin) > 0
      tr = varargin{1};
  else 
      tr = 2000;
  end
  

% prepare structuring element
  st = zeros(3,3,3);
  st(2,:,:) = 1;
  st(:,2,:) = 1;
  st(:,:,2) = 1;
 
% compute the features for each independent component
  total_ics = 0;
  for j = 1:n_of_timeseries
    ica_directory = icadir{j}; 
    len = length(ica_directory);
    if ica_directory(len) == '/'
      ica_directory = ica_directory;
    else
      ica_directory = strcat(ica_directory,'/');
      len = len + 1;
    end
    disp(['Processing ', ica_directory])
    [img,sig,sigFT,mask] = fmri_readicanii(ica_directory);
    if length(varargin) > 1 
      if ~isempty(varargin{2})
        flimits = varargin{2};
        flimits = flimits*((tr/1000)*length(sig));
      else
        flimits = [0.005 0.01 0.095];
        flimits = flimits*((tr/1000)*length(sig));
      end
    else
      flimits = [0.005 0.01 0.095];
      flimits = flimits*((tr/1000)*length(sig));
    end
    flimits(1) = round(flimits(1) + 1);
    flimits(2) = ceil(flimits(2));
    flimits(3) = round(flimits(3));
    sz = size(img);
    n_of_ics{j} = sz(4);
    total_ics = total_ics + n_of_ics{j};
    mask = double(mask);
    img = double(img);
    bsmask = mask - erode3D(mask,st);
    bsind = find(bsmask(:) > 0);
    brain = erode3D(mask,ones(3,3,3));
    brainind = find(brain(:) > 0);
    for i = 1:sz(4)
      features(total_ics - n_of_ics{j} + i,:) = fmri_features4(img(:,:,:,i),sigFT(i,:),sig(i,:),mask,bsind,brainind,blocked,flimits);
    end
  end % for j = 1:n_of_timeseries
  disp('Generating training data')
% compile the classes for training ICs.
  ic_classes = zeros(total_ics,1);
  total_ics = 0;
  for j = 1:n_of_timeseries
    acc{j} = setdiff(1:n_of_ics{j},rej{j});
    ic_classes(total_ics + rej{j}) = 1;
    back_rej{j} = total_ics + rej{j}; % for generating rejected_clusters matrix 
    total_ics = total_ics + n_of_ics{j}; 
  end  
  class0 = find(ic_classes == 0)'; 
% cluster rejected components into 4 classes if the clustering is not provided
  if length(varargin) < 3
    ic_classes_clustered = fmri_cluster_tsamples(features,ic_classes);
  else
    ic_classes_clustered = zeros(length(ic_classes),1);
    for j = 1:n_of_timeseries
        input_rejected_clusters = varargin{3};
        disp('here')
        ic_classes_clustered(back_rej{j}) = input_rejected_clusters{j}(:,2);
    end
  end
  class1 = find(ic_classes_clustered == 1)';
  class2 = find(ic_classes_clustered == 2)';
  class3 = find(ic_classes_clustered == 3)';
  class4 = find(ic_classes_clustered == 4)';
  disp(['Signal class has ',int2str(length(class0)), ' samples']);
  disp(['Noise class 1 has ',int2str(length(class1)), ' samples']);
  disp(['Noise class 2 has ',int2str(length(class2)), ' samples']);
  disp(['Noise class 3 has ',int2str(length(class3)), ' samples']);
  disp(['Noise class 4 has ',int2str(length(class4)), ' samples']);
% generate rejected_clusters
  for j = 1:n_of_timeseries
    rejected_clusters{j}(:,1) = rej{j};
    rejected_clusters{j}(:,2) = ic_classes_clustered(back_rej{j});
  end

% train the classifier
  disp('Training the classifier. This may take a while...');
  global_decision_tree_classifier = train_fmri_multithr(features,class0,class1,class2,class3,class4,np_thr,np_testlevels)

% *************************************
% Classifier evaluation subfunction
% *************************************

function [classes,totalclasses] = classify_fmri_gdtv2(features,tra,trb,trc,trd)

  sz = size(features);

  % each row (of 4) gives the classification result to accepted (0) 
  % or rejected (1) corresponding individual classifiers 

  totalclasses = zeros(sz(1),4);
  classes = zeros(sz(1),1);

  totalclasses(:,1) = (features(:,2) < tra(1)) & (features(:,4) < tra(2)) & (features(:,6) < tra(3));
  totalclasses(:,2) = (features(:,1) < trb(1)) & ((features(:,2) < trb(2)) | (features(:,3) < trb(3)));
  totalclasses(:,3) = (features(:,2) < trc(1)) & (features(:,4) < trc(2)) & (features(:,6) < trc(3));
  totalclasses(:,4) = (features(:,2) < trd(1)) & (features(:,4) < trd(2)) & (features(:,5) < trd(3));
  
  classes = sum(totalclasses,2);



% Train global decision tree classifier for the fmri independent component classification
% *************************************************************************************
% Input : features   :  a matrix containing the (6-component) feature vector for each sample
%          class0    :  indeces of samples in the 'accepted' class
%          class1    :  indeces of samples in the  rejected class A
%          class2    :  indeces of samples in the  rejected class B
%          class3    :  indeces of samples in the  rejected class C
%          class4    :  indeces of samples in the  rejected class D
%          np_thr    :  Neyman Pearson thresholds for the classifier. 
%                       This function generates classifier for these thresholds
%          np_testlevels (optional) : The Neyman Pearson thresholds for component classifiers    
% Output : the global decision tree classifier(s)
% ************************************************************************************* 

function global_decision_tree_classifier = train_fmri_multithr(features,class0,class1,class2,class3,class4,np_thr,np_testlevels)

  if isempty(np_testlevels)
    np_testlevels = [0.0100 0.0150 0.0200 0.0300 0.0400 0.0500 0.0700 0.0900]
  end

  % Neyman Pearson threshols
  global_decision_tree_classifier.np_thr = np_thr;
  % tested Neyman Pearson levels for component classifiers 
  global_decision_tree_classifier.np_testlevels = np_testlevels;

% first compute the allocated np rate: 

  sz0 = length(class0);
  szc(1) = length(class1);
  szc(2) = length(class2);
  szc(3) = length(class3);
  szc(4) = length(class4);



  for i = 1:length(np_testlevels)
   
    [tra_tl{i} tn(i,1) fa(i,1)] = global_decision_tree_v2(features([class0 class1],[2 4 6]),1:sz0,(sz0 + 1):(sz0 + szc(1)),np_testlevels(i),1,[]);
   
    [trb_tl{i} tn(i,2) fa(i,2)] = global_decision_tree_v2(features([class0 class2],[1 2 3]),1:sz0,(sz0 + 1):(sz0 + szc(2)),np_testlevels(i),3,[2 3 1]);
    [trc_tl{i} tn(i,3) fa(i,3)] = global_decision_tree_v2(features([class0 class3],[2 4 6]),1:sz0,(sz0 + 1):(sz0 + szc(3)),np_testlevels(i),1,[]);
    [trd_tl{i} tn(i,4) fa(i,4)] = global_decision_tree_v2(features([class0 class4],[2 4 5]),1:sz0,(sz0 + 1):(sz0 + szc(4)),np_testlevels(i),1,[]);
  end

  total_fa = zeros(length(np_testlevels),length(np_testlevels),length(np_testlevels),length(np_testlevels));
  total_tn = zeros(length(np_testlevels),length(np_testlevels),length(np_testlevels),length(np_testlevels));
  negatives = length(class1) + length(class2) + length(class3) + length(class4);

  max_tn = zeros(length(np_thr),1);
  max_fa = zeros(length(np_thr),1);
  for m = 1:length(np_thr)
    max_tn_indx{m} = [1 1 1 1];
  end

  for i = 1:length(np_testlevels)
    for j = 1:length(np_testlevels)
    %  disp([i j])
      for k = 1:length(np_testlevels)
        for l = 1:length(np_testlevels)
          [classes,totalclasses] = classify_fmri_gdtv2(features,tra_tl{i},trb_tl{j},trc_tl{k},trd_tl{l}); 
          total_fa(i,j,k,l) = sum(classes(class0) > 0)/length(class0);
          total_tn(i,j,k,l) = (sum(classes(class1) > 0) + sum(classes(class2) > 0) + sum(classes(class3) > 0) + sum(classes(class4) > 0))/negatives;
          for m =1:length(np_thr)
            if total_fa(i,j,k,l) < np_thr(m)
              if total_tn(i,j,k,l) > max_tn(m) 
                max_tn(m) = total_tn(i,j,k,l);
                max_tn_indx{m} = [i j k l]; 
                max_fa(m) = total_fa(i,j,k,l);
              end
            end
          end 
        end
      end
    end
  end

  for m = 1:length(np_thr)
    global_decision_tree_classifier.classifier{m}.tra = tra_tl{max_tn_indx{m}(1)};
    global_decision_tree_classifier.classifier{m}.trb = trb_tl{max_tn_indx{m}(2)};
    global_decision_tree_classifier.classifier{m}.trc = trc_tl{max_tn_indx{m}(3)};
    global_decision_tree_classifier.classifier{m}.trd = trd_tl{max_tn_indx{m}(4)};
  end

  global_decision_tree_classifier.truenegatives = max_tn;
  global_decision_tree_classifier.falsealarm = max_fa;
  global_decision_tree_classifier.component_truenegatives = total_tn;
  global_decision_tree_classifier.component_falsealarm = total_fa;
  for m = 1:length(np_thr)
    if max_tn(m) == 0
      disp(strcat('Warning: the classifier at NP =', num2str(np_thr(m)), 'is not going to be useful. Use a higher NP threshold, add training samples, or adjust parameters for classifier training'))
   end
 end 

% **************************************************************************************
% Reads the component maps, IC time courses, IC power spectra, and the brain mask from Melodic output
% requires : Tools for NIfTI MRI by Jimmy Shen 
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=8797&objectType=FILE
% **************************************************************************************  
  
function [img,sig,sigFT,mask] = fmri_readicanii(directory)

  % The file names. You may need to change these to correspond to
% your settings of the melodic.

ICmapfn = 'melodic_IC';
maskfn = 'mask';
ICtcfn = 'melodic_mix';
ICtcftfn = 'melodic_FTmix';


len = length(directory);
if directory(len) == '/'
  ICdir = directory;
else
  ICdir = strcat(directory,'/');
end

hdr_zipped = 0;
if exist(strcat(ICdir,ICmapfn,'.nii'),'file') 
  niifile = load_untouch_nii(strcat(ICdir,ICmapfn,'.nii'));
elseif exist(strcat(ICdir,ICmapfn,'.nii.gz'),'file')
  unix(['gzip -d ',strcat(ICdir,ICmapfn,'.nii.gz')]);
  niifile = load_untouch_nii(strcat(ICdir,ICmapfn,'.nii'));
  unix(['gzip ',strcat(ICdir,ICmapfn,'.nii')]);
elseif exist(strcat(ICdir,ICmapfn,'.img.gz'),'file')
  unix(['gzip -d ',strcat(ICdir,ICmapfn,'.img.gz')]); 
  if exist(strcat(ICdir,ICmapfn,'.hdr.gz'),'file') 
    hdr_zipped = 1;
    unix(['gzip -d ',strcat(ICdir,ICmapfn,'.hdr.gz')]); 
  end
  niifile = load_untouch_nii(strcat(ICdir,ICmapfn));
  unix(['gzip ',strcat(ICdir,ICmapfn,'.img')]); 
  if hdr_zipped 
    unix(['gzip ',strcat(ICdir,ICmapfn,'.hdr')]); 
  end
else
  niifile = load_untouch_nii(strcat(ICdir,ICmapfn));
end
img = niifile.img;

hdr_zipped = 0;
if exist(strcat(ICdir,maskfn,'.nii'),'file') 
  niifile = load_untouch_nii(strcat(ICdir,maskfn,'.nii'));
elseif exist(strcat(ICdir,maskfn,'.nii.gz'),'file')
  unix(['gzip -d ',strcat(ICdir,maskfn,'.nii.gz')]);
  niifile = load_untouch_nii(strcat(ICdir,maskfn,'.nii'));
  unix(['gzip ',strcat(ICdir,maskfn,'.nii')]);
elseif exist(strcat(ICdir,maskfn,'.img.gz'),'file')
  unix(['gzip -d ',strcat(ICdir,maskfn,'.img.gz')]); 
  if exist(strcat(ICdir,maskfn,'.hdr.gz'),'file') 
    hdr_zipped = 1;
    unix(['gzip -d ',strcat(ICdir,maskfn,'.hdr.gz')]); 
  end
  niifile = load_untouch_nii(strcat(ICdir,maskfn));
  unix(['gzip ',strcat(ICdir,maskfn,'.img')]); 
  if hdr_zipped 
    unix(['gzip ',strcat(ICdir,maskfn,'.hdr')]); 
  end
else
  niifile = load_untouch_nii(strcat(ICdir,maskfn));
end
mask = niifile.img;
sz = size(img);

fid = fopen(strcat(ICdir,ICtcfn),'rt');
[sig,count] = fscanf(fid,'%f');
sig = reshape(sig,sz(4),count/sz(4));
fclose(fid);
fid = fopen(strcat(ICdir,ICtcftfn),'rt');
[sigFT,count] = fscanf(fid,'%f');
sigFT = reshape(sigFT,sz(4),count/sz(4));
fclose(fid);

% Clustering the rejected ICs in the training data 

function [class,score,class_score] = fmri_cluster_tsamples(features,status)

  sz = size(status);
  class = zeros(sz(1),1);
  score = zeros(sz(1),6);
  class_score = zeros(sz(1),4);

  acc_blocked = find(status(:) == 0);

  acblen = length(acc_blocked);


  for i = 1:length(status)
    if status(i) == 1
      for j = 1:6
        score(i,j) = sum(features(acc_blocked,j) < features(i,j));
      end

      if (score(i,1) < acblen*0.025) | (((score(i,4) - score(i,1))/(score(i,4) + 1)) > 0.5)
        class(i) = 2;
      elseif ((score(i,2) - score(i,3))/(score(i,2) + 1) > 0.75) & ...
	((score(i,4) - score(i,1))/(score(i,4) + 1) > 0)
        class(i) = 2; 
      elseif (score(i,5) < acblen*0.025) | ((score(i,6) - score(i,5))/(score(i,6) + 1) > 0.75)
        class(i) = 4;
      elseif (score(i,5) - score(i,6))/(score(i,5) + 1) > 0.75
        class(i) = 3;
      else 
        class(i) = 1;
      end
    end
  end
  

