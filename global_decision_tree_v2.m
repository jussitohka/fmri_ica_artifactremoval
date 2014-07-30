% Finds the global decision tree decision boundary in the Neyman Pearson setting
% Copyright (C) 2005 - 2007 Jussi Tohka, Institute of Signal Processing, Tampere University of
% Technology, P.O. Box 553, FIN-33101, Finland
% E-mail: jussi.tohka@tut.fi
%
% The method is described in:
% J. Tohka, K. Foerde, A.R. Aron, S. M. Tom, A.W. Toga, and
% R.A. Poldrack. 
% Automatic Independent Component Labeling for Artifact Removal in fMRI.
% NeuroImage, In press, 2007. http://dx.doi.org/10.1016/j.neuroimage.2007.10.013
% ----------------------------------------------------------------------------
% Permission to use, copy, modify, and distribute this software 
% for any purpose and without fee is hereby
% granted, provided that the above copyright notice appear in all
% copies.  The author and Tampere University of Technology make no representations
% about the suitability of this software for any purpose.  It is
% provided "as is" without express or implied warranty.
% ----------------------------------------------------------------------------
%  
% 
%  decision rules : (x1 < tr1 AND/OR x2 < tr2) AND/OR x3 < tr3
%                 :  x1 < tr1 AND/OR (x2 < tr2 AND/OR x3 < tr3)
%                 : (x1 < tr1 AND/OR x3 < tr3) AND/OR x2 < tr2

% scale features so that the features in class 1 are greater than in the class 2
% nplimit gives the percentage of the class1 features to be classified to the class 2
% decrule gives the decision rule
% decrule = 1 : AND AND
% decrule = 2 : AND OR
% decrule = 3 : OR AND
% decrule = 4 : OR OR
% with feat order the features can be re-ordered
% decision rule is always of the form  (x1 < tr1 AND/OR x2 < tr2) AND/OR x3 < tr3 
% if this is not required give the empty matrix

function [tr,tnrate,farate] = global_decision_tree_v2(features,class1,class2,nplimit,decrule,feat_order)

if ~isempty(feat_order)
  features = features(:,feat_order);
end

sz = size(features);
sz1 = length(class1);
sz2 = length(class2);

nptr = max(floor(nplimit*sz1),2);



% it is assumed that sz2 < sz1
for i = 1:sz(2)
  sorted_features2(:,i) = sort(features(class2,i));
  sorted_features1(:,i) = sort(features(class1,i));
  sorted_features(:,i) = sort(features(:,i));
end
features1 = features(class1,:);
features2 = features(class2,:); 

if decrule == 1
  up_testtr = sorted_features1(sz1,:);
  low_testtr = sorted_features1(nptr - 1,:) - 0.001;
end
if decrule == 2
  up_testtr = sorted_features1(sz1,:);
  low_testtr = sorted_features1(1,:) - 0.001;
end
if decrule == 3
  up_testtr = sorted_features1(sz1,:);
  low_testtr = sorted_features1(1,:) - 0.001;
  low_testtr(3) = sorted_features1(nptr - 1,3) - 0.001;
end
if decrule == 4
  up_testtr = sorted_features1(nptr + 2,:);
  low_testtr = sorted_features1(1,:) - 0.001;
end



for i = 1:sz(2)
  testind = find((sorted_features2(:,i) >= low_testtr(i)) & (sorted_features2(:,i) <= up_testtr(i)));
  if length(testind) > 1 
    test_tr{i} = sorted_features2(testind,i) + 0.00001;
    testnumber(i) = length(testind);
  else
    test_tr{i} = [sorted_features2(testind,i) low_testtr(i)] + 0.00001;
    testnumber(i) = length(testind) + 1;
  end
end


tr = [1 1 1];
fa_max = 0;
tn_max = 0;

switch decrule
  case 1  % AND AND
    for i = 1:testnumber(1)
      
      fa = length(find((features1(:,1) < test_tr{1}(i)) & (features1(:,2) < test_tr{2}(1)) ...
                                                                          & (features1(:,3) < test_tr{3}(1))));
      tn = length(find((features2(:,1) < test_tr{1}(i)) & (features2(:,2) < test_tr{2}(1)) & ...
                                                                            (features2(:,3) < test_tr{3}(1))));  
      if fa > (nptr + 1)   
        break
      else
        if tn > tn_max 
          fa_max = fa;
          tn_max = tn;
          tr = [i 1 1];
        end  
      end
      for j = 2:testnumber(2)
        fa = length(find((features1(:,1) < test_tr{1}(i)) & (features1(:,2) < test_tr{2}(j)) ...
                                                                          & (features1(:,3) < test_tr{3}(1))));
        tn = length(find((features2(:,1) < test_tr{1}(i)) & (features2(:,2) < test_tr{2}(j)) & ...
                                                                            (features2(:,3) < test_tr{3}(1))));     
        if fa > (nptr + 1)   
          break
        else
          if tn > tn_max 
            fa_max = fa;
            tn_max = tn;
            tr = [i j 1];
          end
        end
        for k = 2:testnumber(3)
          fa = length(find((features1(:,1) < test_tr{1}(i)) & (features1(:,2) < test_tr{2}(j)) ...
                                                                          & (features1(:,3) < test_tr{3}(k))));
          tn = length(find((features2(:,1) < test_tr{1}(i)) & (features2(:,2) < test_tr{2}(j)) & ...
                                                                            (features2(:,3) < test_tr{3}(k))));
          if fa > (nptr + 1)   
            break;
          else
            if tn > tn_max 
              fa_max = fa;
              tn_max = tn;
              tr = [i j k];
            end

          end
        end
      end
    end
  case 2  % AND OR
    for i = 1:testnumber(1)
      
      fa = length(find( ((features1(:,1) < test_tr{1}(i)) & (features1(:,2) < test_tr{2}(1))) ...
                                                                          | (features1(:,3) < test_tr{3}(1))));
      tn = length(find( ((features2(:,1) < test_tr{1}(i)) & (features2(:,2) < test_tr{2}(1))) | ...
                                                                            (features2(:,3) < test_tr{3}(1))));  
      if fa > (nptr + 1)   
        break
      else
        if tn > tn_max 
          fa_max = fa;
          tn_max = tn;
          tr = [i 1 1];
        end 
      end
      for j = 2:testnumber(2)
        fa = length(find( ((features1(:,1) < test_tr{1}(i)) & (features1(:,2) < test_tr{2}(j))) ...
                                                                          | (features1(:,3) < test_tr{3}(1))));
        tn = length(find( ((features2(:,1) < test_tr{1}(i)) & (features2(:,2) < test_tr{2}(j))) | ...
                                                                            (features2(:,3) < test_tr{3}(1))));     
        if fa > (nptr + 1)   
          break
        else
          if tn > tn_max 
            fa_max = fa;
            tn_max = tn;
            tr = [i j 1];
          end 
        end
        for k = 2:testnumber(3)
          fa = length(find( ((features1(:,1) < test_tr{1}(i)) & (features1(:,2) < test_tr{2}(j))) ...
                                                                          | (features1(:,3) < test_tr{3}(k))));
          tn = length(find( ((features2(:,1) < test_tr{1}(i)) & (features2(:,2) < test_tr{2}(j))) | ...
                                                                            (features2(:,3) < test_tr{3}(k))));
          if fa > (nptr + 1)   
            break;
          else
            if tn > tn_max 
              fa_max = fa;
              tn_max = tn;
              tr = [i j k];
            end 
          end
        end
      end
    end     
  case 3  % OR AND
    for i = 1:testnumber(1)
      
      fa = length(find( ((features1(:,1) < test_tr{1}(i)) | (features1(:,2) < test_tr{2}(1))) ...
                                                                          & (features1(:,3) < test_tr{3}(1))));
      tn = length(find( ((features2(:,1) < test_tr{1}(i)) | (features2(:,2) < test_tr{2}(1))) & ...
                                                                            (features2(:,3) < test_tr{3}(1))));  
      if fa > (nptr + 1)   
        break
      else
        if tn > tn_max 
          fa_max = fa;
          tn_max = tn;
          tr = [i 1 1];
        end 
      end
      for j = 2:testnumber(2)
        fa = length(find( ((features1(:,1) < test_tr{1}(i)) | (features1(:,2) < test_tr{2}(j))) ...
                                                                          & (features1(:,3) < test_tr{3}(1))));
        tn = length(find( ((features2(:,1) < test_tr{1}(i)) | (features2(:,2) < test_tr{2}(j))) & ...
                                                                            (features2(:,3) < test_tr{3}(1))));     
        if fa > (nptr + 1)   
          break
        else
          if tn > tn_max 
            fa_max = fa;
            tn_max = tn;
            tr = [i j 1];
          end 
        end
        for k = 2:testnumber(3)
          fa = length(find( ((features1(:,1) < test_tr{1}(i)) | (features1(:,2) < test_tr{2}(j))) ...
                                                                          & (features1(:,3) < test_tr{3}(k))));
          tn = length(find( ((features2(:,1) < test_tr{1}(i)) | (features2(:,2) < test_tr{2}(j))) & ...
                                                                            (features2(:,3) < test_tr{3}(k))));
          if fa > (nptr + 1)   
            break;
          else
            if tn > tn_max 
              fa_max = fa;
              tn_max = tn;
              tr = [i j k];
            end 
          end
        end
      end
    end     
  case 4  % OR OR
    for i = 1:testnumber(1)
      
      fa = length(find( ((features1(:,1) < test_tr{1}(i)) | (features1(:,2) < test_tr{2}(1))) ...
                                                                          | (features1(:,3) < test_tr{3}(1))));
      tn = length(find( ((features2(:,1) < test_tr{1}(i)) | (features2(:,2) < test_tr{2}(1))) | ...
                                                                            (features2(:,3) < test_tr{3}(1))));  
      if fa > (nptr + 1)   
        break
      else
        if tn > tn_max 
          fa_max = fa;
          tn_max = tn;
          tr = [i 1 1];
        end 
      end
      for j = 2:testnumber(2)
        fa = length(find( ((features1(:,1) < test_tr{1}(i)) | (features1(:,2) < test_tr{2}(j))) ...
                                                                          | (features1(:,3) < test_tr{3}(1))));
        tn = length(find( ((features2(:,1) < test_tr{1}(i)) | (features2(:,2) < test_tr{2}(j))) | ...
                                                                            (features2(:,3) < test_tr{3}(1))));     
        if fa > (nptr + 1)   
          break
        else
          if tn > tn_max 
            fa_max = fa;
            tn_max = tn;
            tr = [i j 1];
          end 
        end
        for k = 2:testnumber(3)
          fa = length(find( ((features1(:,1) < test_tr{1}(i)) | (features1(:,2) < test_tr{2}(j))) ...
                                                                          | (features1(:,3) < test_tr{3}(k))));
          tn = length(find( ((features2(:,1) < test_tr{1}(i)) | (features2(:,2) < test_tr{2}(j))) | ...
                                                                            (features2(:,3) < test_tr{3}(k))));
          if fa > (nptr + 1)   
            break;
          else
            if tn > tn_max 
              fa_max = fa;
              tn_max = tn;
              tr = [i j k];
            end 
          end
        end
      end
    end
  end

if sz2 > 0
  tnrate = tn_max/sz2;
else
  tnrate = tn_max;
end
farate = fa_max/sz1;

for i = 1:sz(2)
  tr(i) = test_tr{i}(tr(i));
end
 
if ~isempty(feat_order)
   tr(feat_order) = tr;
end


