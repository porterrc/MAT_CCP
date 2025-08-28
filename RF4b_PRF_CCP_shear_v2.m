 close all; clear all;
addpath ./Functions
addpath ./Data
addpath ./Functions/Color_Palette

decon = 'ID' % Choose Deconvolution method

lat1 = 34.55 % southern lat for grid
lat2 = 36.1 % northern lat for grid
lon1 = -112.6 % western edge for grid 
lon2 = -110.7 % eastern edge for grid
dzi = 0.5 % depth increment 
z_max = 100 % max depth
bin_sz = 2 % bin size in km  
bin_sm = 3 %1number of bins to smooth from % og value =5
gauss_width = 6; %std in km for gauss weighting using normpdf % og value 20
shear_file = 'WUS324.r0.0.nc'; %'FWT-SouthAmerica-2022.r0.0.nc' 
vmodel = 'Wolf.tvel'
minv = 2.5; %minimum velocity to accept in model, really low values produce unrealistic results
max_amp = 2; % remove traces with amps great than this value
min_amp = -2; % remove traces with amps less than this value
blim = [0 360]
boot_strap = 20 % number of bootstraps to calculate
xdepth = [22 24];
xcenter = [35.33 -111.69];
xradius = [16];
bin_min = 4; 
bin_max = 15; 
min_traces = 25;

load RF_5_Merge.mat

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

fprintf('Converting Shear Model to km grid \n')
[sx sy sz svmodel pvmodel] = shear_loc_FWT(shear_file,lat1,lat2,lon1,lon2,dd);

fprintf('Tracing Rays and Locate on Grid \n')
%[z, xpierce, ypierce, seis] = ray_trace_prf_shear_v3b(Data,sname,stnx,stny,stel,dzi,z_max,sx,sy,sz,svmodel,vmodel,minv,max_amp,min_amp,blim);
[z, xpierce, ypierce, baz,seis,tnames] = ray_trace_prf_ps(Data,sname,stnx,stny,stel,dzi,z_max,sx,sy,sz,svmodel,pvmodel,vmodel,minv,max_amp,min_amp,blim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Convert X and Ys to Lat Lon \n');
lat = lat1+km2deg(y);
lon = lon1+x./dd; 
xlpierce = lon1 + xpierce./dd;
ylpierce = lat1 + km2deg(ypierce,'earth');
[lats lons] = meshgrid(lat,lon);

fprintf('Removing Traces in Exclusion Zone \n')
%[xpierce, ypierce,ypierce,xlpierce, seis] = exclude_traces(z, xpierce, ypierce,xlpierce,ylpierce, seis,xcenter,xdepth,xradius)

    % ind = find( z>xdepth(1) & z< xdepth(2));
    % for i = 1:length(ind)
    %     [~,d] = knnsearch(xcenter,[ylpierce(ind(i),:)',xlpierce(ind(i),:)']);
    %     dind = find(d*deg2km(1) < xradius);
    %     xpierce(ind(i):end,dind) = nan;
    %     ypierce(ind(i):end,dind) = nan;
    %     xlpierce(ind(i):end,dind) = nan;
    %     ylpierce(ind(i):end,dind) = nan;
    %     seis(ind(i):end,dind) = nan;
    % end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Plotting Stations, Velocities, Gridpoints, and Piercing Points at 10 km depth\n')
plot_array(X,Y,z,stnx,stny,sname,xpierce,ypierce,sx,sy,sz,svmodel,10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Create Data Grid\n')
% [CCP CCP_RAW BAZ_RAW HITS BSTD] = Grid_data_anisotropy(xpierce,ypierce,seis,baz,x,y,z,bin_sz,bin_sm,gauss_width,boot_strap);
% [CCPiso CCPan CCPaz] = CCP_anistroy(CCP,CCP_RAW,BAZ_RAW);

%[CCP HITS BSTD] = Grid_data_v3(xpierce,ypierce,seis,x,y,z,bin_sz,bin_sm,gauss_width,boot_strap);
[CCP HITS BSTD GSIZE] = Grid_data_Variable(xpierce,ypierce,seis,x,y,z,bin_min,bin_max,min_traces,gauss_width,boot_strap,tnames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Convert X and Ys to Lat Lon \n');
lat = lat1+km2deg(y);
lon = lon1+x./dd; 
xlpierce = lon1 + xpierce./dd;
ylpierce = lat1 + km2deg(ypierce,'earth');
[lats lons] = meshgrid(lat,lon);

if exist('GSIZE')
    save('CCP_5s.mat','CCP','lat','lon','x','y','z','Stns','xpierce','ypierce','xlpierce','ylpierce','HITS','BSTD','GSIZE','bin_sz','bin_sm')
else
    save('CCP_5s.mat','CCP','lat','lon','x','y','z','Stns','xpierce','ypierce','xlpierce','ylpierce','HITS','BSTD','bin_sz','bin_sm')
end
save('CCP_ray_info_5s.mat','lat','lon','x','y','z','Stns','xpierce','ypierce','xlpierce','ylpierce','seis','Stns')
%save('CCP_5_Anis.mat','CCP','CCPiso','CCPan','CCPaz','BAZ_RAW','CCP_RAW','lat','lon','x','y','z','Stns','xpierce','ypierce','xlpierce','ylpierce','HITS','BSTD','bin_sz','bin_sm')
