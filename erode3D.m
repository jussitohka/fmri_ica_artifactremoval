% Copyright (C) 2001 Jussi Tohka, Institute of Signal Processing, Tampere University of
% Technology, P.O. Box 553, FIN-33101, Finland
% E-mail: jussi.tohka@tut.fi
% ----------------------------------------------------------------------------
% Permission to use, copy, modify, and distribute this software 
% for any purpose and without fee is hereby
% granted, provided that the above copyright notice appear in all
% copies.  The author and Tampere University of Technology make no representations
% about the suitability of this software for any purpose.  It is
% provided "as is" without express or implied warranty.
% ----------------------------------------------------------------------------
% Three dimensional morphological (binary) erosion.
% Input: X is the image to be processed. X(i,j,k) == 1 if voxel
% (i,j,k) \in X otherwise X(i,j,k) == 0.
% B is the structuring element. Note that the size of B should be o1 x 
% o2 x o3 where o1, o2, o3 are odd numbers. The origin is (using previous 
% notation) B(ceil(o1),ceil(o2),ceil(o3)) i.e the centre voxel.

% Note: Output XdB of the same size than X. If this is not what you 
% are after then change the call to convn.

function XeB = erode3D(X,B)
  s = size(B);
  if rem(s(1),2) == 0 | rem(s(2),2) == 0 | rem(s(3),2) == 0
    disp('Improper size of structuring element. Doublecheck the output.');
  end
  Belements = sum(sum(sum(B)));   
  XeB = convn(X,B,'same');
  XeB = XeB == Belements;
  