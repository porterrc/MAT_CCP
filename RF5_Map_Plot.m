close all; clear all;
addpath ./Functions
addpath ./Data
addpath ./Functions/Color_Palette
%
%depth = 9; %depth below sea level in km
search_polarity =-1; %postive if looking for max, negative if looking for min

 depth_min = 5; %depth below sea level in km
 depth_max = 18; %depth below sea level in km
  %depth_min = 30; %depth below sea level in km
  %depth_max = 58; %depth below sea level in km


load CCP_5s.mat
vents = load('SFVF_vents.txt');

%ind = find(abs(z-depth) == min(abs(z-depth)));
indmin = find(abs(z-depth_min) == min(abs(z-depth_min)));
indmax = find(abs(z-depth_max) == min(abs(z-depth_max)));

[Lon,Lat] = meshgrid(lon,lat);
%RF_Amp = -squeeze(CCP(ind,:,:));
RF_Temp = -1.*search_polarity*squeeze(CCP(indmin:indmax,:,:));
[RF_Amp RF_Ind] = ((max(RF_Temp,[],1)));
RF_Amp = squeeze(RF_Amp.*search_polarity); RF_Ind = squeeze(RF_Ind);
IDepth = (z(RF_Ind) + depth_min);

[pks,ws,ps,dpk] = Moho_Picker(-1.*search_polarity*CCP,lon,lat,z,depth_min,depth_max);
pk_amp = search_polarity.*squeeze(max(pks,[],1));
pk_amp(pk_amp == 0) = nan;

s_list = fields(Stns);
for i = 1:length(s_list)
    slon(i) = Stns.(char(s_list(i))).Station_Data.Longitude;
    slat(i) = Stns.(char(s_list(i))).Station_Data.Latitude;
    selv(i) = Stns.(char(s_list(i))).Station_Data.Elevation;
end

sboundary = boundary(slon',slat');
sta_in = inpolygon(Lon,Lat,slon(sboundary)',slat(sboundary)');

pk_amp(~sta_in) = nan;
dpk(~sta_in) = nan;

latq = [min(lat):.01:max(lat)];
lonq = [min(lon):.01:max(lon)];
[Lonq,Latq] = meshgrid(lonq,latq);
RF_Ampq = interp2(Lon,Lat,pk_amp,Lonq,Latq);
Idepthq = interp2(Lon,Lat,dpk,Lonq,Latq);

%[SA,mlon,mlat] = TopoDownload([min(lon),max(lon)],[min(lat),max(lat)],'default');





 Stk = ncread('Pcb.nc','PcB');
 Slon = ncread('Pcb.nc','lon');
 Slat = ncread('Pcb.nc','lat'); 
 Stk(Stk == -9999) = nan;
 
 [Sy Sx] = meshgrid(Slat,Slon);

 St = interp2(Sy,Sx,Stk,Latq,Lonq,'linear',nan);
%%%%%%%%%%%Set Color Map%%%%%%%%
       load No_Green
       cprf = no_green;
       %colormap(cprf); sca = 0.5; caxis([-sca sca]); 
%%%%%%%%%

% figure
% ax1 = axes;
% %worldmap([min(lat),max(lat)],[min(lon),max(lon)])
%  M = surf(ax1,mlon,mlat,SA)
%  set(M,'EdgeColor','none')
%  MCmap = demcmap(SA)
%  view(2)
% caxis([2.5 4])
%surfm(dlat,dlon,svmodel001(:,:,30))
%load coastlines 
%geoshow(coastlat,coastl
%worldmap([min(lat),max(lat)],[min(lon),max(lon)])
%%Hide the top axes

% ax2 = axes;
% hold on
% C = surf(ax2,Lon,Lat,IDepth,'FaceAlpha',0.75)
% set(C,'EdgeColor','none')
% %caxis([-1 1])
% view(2)
% 
% linkaxes([ax1,ax2])
% ax2.Visible = 'off';
% ax2.XTick = [];
% ax2.YTick = [];
% 
% colormap(ax1,gray)
% %colormap(ax2,cprf)
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);
% %colormap(flipud(colormap))
% colorbar
% %caxis([-1 1])
% 
% xlim([min(lon),max(lon)])
% ylim([min(lat),max(lat)])
% 
% s_list = fields(Stns);
% for i = 1:length(s_list)
%     scatter(ax2,Stns.(char(s_list(i))).Station_Data.Longitude,Stns.(char(s_list(i))).Station_Data.Latitude,5,'red','filled');
% end
% 
% figure
% worldmap([min(lat),max(lat)],[min(lon),max(lon)])
% C = surfm(Lat,Lon,RF_Amp)
% set(C,'EdgeColor','none')
% caxis([-.5 .5])
% colormap(cprf)
% colorbar
% view(2)
cones = [35.342663, -112.006069; 35.408056, -111.850833; 35.200133, -112.205213; 35.333333, -111.666667; 34.75, -111.6; 34.96 -111.5; 34.8 -111.85; 34.53 -111.73];



figure
worldmap([min(lat),max(lat)],[min(lon),max(lon)])
hold on
C = surfm(Latq,Lonq,Idepthq)
set(C,'EdgeColor','none')
caxis([depth_min depth_max])
scatterm(vents(:,2),vents(:,1),15,'filled','k')
scatterm(cones(:,1),cones(:,2),40,'filled','r')
%colormap(cprf)
colorbar
view(2)

figure
worldmap([min(lat),max(lat)],[min(lon),max(lon)])
hold on
load('Sunset_Swarm.mat')
C = surfm(Latq,Lonq,RF_Ampq)
set(C,'EdgeColor','none')
scatterm(vents(:,2),vents(:,1),15,'filled','k')
scatterm(cones(:,1),cones(:,2),40,'filled','r')
%scatterm(Event_Data.PreferredLatitude,Event_Data.PreferredLongitude,'filled')
%caxis([0.1 1])
caxis([-.75 .75])
colormap(cprf)
colorbar
view(2)

figure
hold on
[GLon,GLat,Gvals] = Grav_Data('SFVF_Grav.nc');
worldmap([min(lat),max(lat)],[min(lon),max(lon)])
C = surfm(GLat,GLon,Gvals)
set(C,'EdgeColor','none')
%caxis([0.1 1])
caxis([-250 -130])
colormap(cprf)
scatterm(vents(:,2),vents(:,1),'filled','k')
scatterm(cones(:,1),cones(:,2),40,'filled','r')
colorbar
view(2)


figure
hold on
worldmap([min(lat),max(lat)],[min(lon),max(lon)])
C = surfm(Latq,Lonq,St*0.3048/1000)
set(C,'EdgeColor','none')
%caxis([0.1 1])
%caxis([-250 -130])
colormap(cprf)
scatterm(vents(:,2),vents(:,1),'filled','k')
scatterm(cones(:,1),cones(:,2),40,'filled','r')
colorbar
view(2)


figure
hold on
worldmap([min(lat),max(lat)],[min(lon),max(lon)])
[blat, blon, Bt] = Basement_Depth_Map(min(lat),max(lat),min(lon),max(lon));
C = surfm(blat,blon,smoothdata2(Bt,'movmean',61))
set(C,'EdgeColor','none')
%caxis([0.1 1])
%caxis([-250 -130])
scatterm(vents(:,2),vents(:,1),'filled','k')
scatterm(cones(:,1),cones(:,2),40,'filled','r')
colorbar
view(2)


figure
worldmap([min(lat),max(lat)],[min(lon),max(lon)])
bdepthq = interp2(Lon,Lat,dpk,blon,blat);
Qt = interp2(blon',blat',Bt',slon,slat);
scatterm(slat,slon,40,Qt,'filled')
brat = bdepthq./Bt;
brat((brat > 20)) = nan;
%C = surfm(blat,blon,brat)
set(C,'EdgeColor','none')
%caxis([0.1 1])
%caxis([-250 -130])
scatterm(vents(:,2),vents(:,1),'filled','k')
scatterm(cones(:,1),cones(:,2),40,'filled','r')
colorbar
view(2)
% grdwrite2(lon,lat,idepth,'RF_depth.grd')
% grdwrite2(lon,lat,iamp,'RF_amp.grd')
% 
 lat_export = reshape(Lat,[],1);
 lon_export = reshape(Lon,[],1);
 rfd_export = reshape(dpk,[],1);
 rfa_export = reshape(pk_amp,[],1);

 fileID = fopen(['RF_Amp_Depth.txt'],'w');
 fprintf(fileID,'%6.2f %6.2f %6.2f %6.2f\n',[lon_export'; lat_export'; rfd_export'; rfa_export']);
 fclose(fileID);

%save('Moho.mat','Lat','Lon','idepth','iamp')