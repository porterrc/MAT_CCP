function [topokm] = topography_angle_v4(x_val,y_val,rf_dist,resolution,dplot);
%[SA,refvec] = gtopo30('~/Dropbox/Data/Topo_Files/',2,[min(y_val)-.5,max(y_val)+.5],[min(x_val)-.5,max(x_val)+.5]);
[SA,lon,lat] = TopoDownload([min(x_val)-.5,max(x_val)+.5],[min(y_val)-.5,max(y_val)+.5],resolution);

i = find(isnan(SA));
SA(i)=zeros(size(i));

sz_SA = size(SA);
X = lon;
Y = lat;

topo = interp2(X',Y',SA',x_val,y_val);
a = size(topo);
topokm = topo/1000;

if ~exist('dplot','var')
 topokm(a+1) = 0;
 topokm(a+2) = 0;
 top_r = rf_dist(1,:);
 top_r(a+1) = top_r(a);
 top_r(a+2) = top_r(1);
 fill(top_r,topokm,[139/255 137/255 137/255])
end
