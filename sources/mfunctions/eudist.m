function y=eudist(loc1,loc2)
%% calculates cord distance bw two sets of points given (x1,y1)(n1x2) and
%% (x2,y2)(n2x2)(in km, etc.)

  if nargin <2 
    loc2 = loc1;
  end  
  size1 = size(loc1);
  size2 = size(loc2);
  n1 = size(loc1,1);
  n2 = size(loc2,1);
  if size1(1,2)==2 & size2(1,2)==2
    x1 = reshape(loc1(:,1),n1,1);
    x2 = reshape(loc2(:,1),n2,1);
    y1 = reshape(loc1(:,2),n1,1);
    y2 = reshape(loc2(:,2),n2,1);
    y = (repmat(x1,1,n2) - repmat(x2',n1,1)).^2 +  ...
        (repmat(y1,1,n2) - repmat(y2',n1,1)).^2;
    y = sqrt(y);
  else
    disp 'loc1 should have 2 columns'
  end
  