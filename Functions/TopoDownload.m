function [z lon lat] = TopoDownload(lonlim, latlim,resolution) 
%resolution can be default, med, high, max

url = ['https://www.gmrt.org:443/services/GridServer?minlongitude=' num2str(lonlim(1)) '&maxlongitude=' num2str(lonlim(2)) '&minlatitude=' num2str(latlim(1)) '&maxlatitude=' num2str(latlim(2)) '&format=netcdf&resolution=' resolution '&layer=topo'];
options = weboptions;
options.Timeout = 60;
websave('topo.nc',url,options);
z = ncread('topo.nc','z');
x_rng = ncread('topo.nc','x_range');
y_rng = ncread('topo.nc','y_range');
spc = ncread('topo.nc','spacing');
sz = ncread('topo.nc','dimension');
lon = x_rng(1):spc(1):x_rng(2);
lat = y_rng(1):spc(2):y_rng(2);
[lat lon] = meshgrid(lat, lon);
z = fliplr(reshape(z,size(lon)));