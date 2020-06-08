function y=eucldist(loc1,loc2)
%% calculates cord distance bw two sets of points given x and y
%% loc1=[x1;y1] (n1x2), loc2=[x2,y2] (n2x2)
  if nargin==1
    loc2 = loc1;
  end
  x1 = loc1(:,1);x2 = loc2(:,1);
  y1 = loc1(:,2);y2 = loc2(:,2);
  n1 = length(x1); n2 = length(x2);
  y = (repmat(x1,1,n2)-repmat(x2',n1,1)).^2 + (repmat(y1,1,n2)- ...
                                               repmat(y2',n1,1)).^2 ;
  y = sqrt(y)*6400*pi/180;