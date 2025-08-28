close all; clear all;
addpath ~/Dropbox/Data/Topo_Files
addpath ./Functions
addpath ./Data
addpath ./Functions/Color_Palette

decon = 'ID' % Choose Deconvolution method

lat1 = 34.55 % southern lat for grid
lat2 = 36.1 % northern lat for grid
lon1 = -112.6 % western edge for grid 
lon2 = -110.7 % eastern edge for grid
% lat1 = 35.1 % southern lat for grid
% lat2 = 35.6 % northern lat for grid
% lon1 = -111.9 % western edge for grid 
% lon2 = -111.15 % eastern edge for grid
dzi = 0.25 % depth increment 
z_max = 100 % max depth
bin_sz = 2 % bin size in km  
bin_sm = 3 % number of bins to smooth from
gauss_width = 6; %std in km for gauss weighting using normpdf
vmodel = 'Wolf.tvel'
max_amp = 1.5; % remove traces with amps great than this value
min_amp = -1.5; % remove traces with amps less than this value
blim = [-360 360]
boot_strap = 50 % number of bootstraps to calculate
bin_min = 2; 
bin_max = 15; 
min_traces = 25;


load RF_5_mat.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Building Grid \n')
[ylen] = deg2km(distance(lat1,lon1,lat2,lon1));
dd = deg2km(mean([distance(lat1,lon1,lat1,lon2) distance(lat2,lon1,lat2,lon2)])/(lon2-lon1)); %km2deg for lat 
[xlen] = dd*(lon2-lon1);
x = 0:bin_sz:xlen+bin_sz; % x values
y = 0:bin_sz:ylen+bin_sz; % y values
[X,Y] = meshgrid(x,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Locating Stations on Grid \n')
 [sname stnx stny stel] = stn_loc_v3(Stns,lat1,lat2,lon1,lon2,dd);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Tracing Rays and Locate on Grid \n')
[z, xpierce, ypierce,baz, seis] = ray_trace_prf_v4(Data,sname,stnx,stny,stel,dzi,z_max,vmodel,max_amp,min_amp,blim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Plotting Stations, Gridpoints, and Piercing Points at 20 km depth\n')
figure
hold on;
plot(X,Y,'+k');
plot(stnx,stny,'*');
text(stnx,stny,sname);
ind = find(z == 20);
%plot3(xpierce(ind,:),ypierce(ind,:),ones(size(xpierce,2))*40,'r.');
%plot3(xpierce,ypierce,z);
set(gca,'Zdir','reverse')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Create Data Grid\n')
%[CCP HITS BSTD] = Grid_data_v3(xpierce,ypierce,seis,x,y,z,bin_sz,bin_sm,gauss_width,boot_strap);
[CCP HITS BSTD GSIZE] = Grid_data_Variable(xpierce,ypierce,seis,x,y,z,bin_min,bin_max,min_traces,gauss_width,boot_strap);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Convert X and Ys to Lat Lon \n');
lat = lat1+km2deg(y);
lon = lon1+x./dd; 
xlpierce = lon1 + xpierce./dd;
ylpierce = lat1 + km2deg(ypierce,'earth');
[lats lons] = meshgrid(lat,lon);


if exist('GSIZE')
    save('CCP_5node.mat','CCP','lat','lon','x','y','z','Stns','xpierce','ypierce','xlpierce','ylpierce','HITS','BSTD','GSIZE','bin_sz','bin_sm')
else
    save('CCP_5.mat','CCP','lat','lon','x','y','z','Stns','xpierce','ypierce','xlpierce','ylpierce','HITS','BSTD','bin_sz','bin_sm')
end
save('CCP_ray_info_5node.mat','lat','lon','x','y','z','Stns','xpierce','ypierce','xlpierce','ylpierce','seis','Stns')