function y=cordist(loc1,loc2)
%% calculates cord distance bw two sets of points given lon and lat
%% loc=[lon;lat], dist is in km

R = 6378.388;
if nargin==1
    loc2 = loc1;
end
coslat1 = cos((loc1(:, 2) * pi)/180);
sinlat1 = sin((loc1(:, 2) * pi)/180);
coslon1 = cos((loc1(:, 1) * pi)/180);
sinlon1 = sin((loc1(:, 1) * pi)/180);
coslat2 = cos((loc2(:, 2) * pi)/180);
sinlat2 = sin((loc2(:, 2) * pi)/180);
coslon2 = cos((loc2(:, 1) * pi)/180);
sinlon2 = sin((loc2(:, 1) * pi)/180);

pp = [coslat1 .* coslon1, coslat1 .* sinlon1, sinlat1] * ...
     [coslat2 .* coslon2, coslat2 .* sinlon2, sinlat2]';
pp(pp>1)=1;

y = R * sqrt( 2* ( 1 - pp ) );