close all; clear all;
addpath ./Functions
addpath ./Data
addpath ./Functions/Color_Palette
%
input_swtch = 1; % Pick x_section location on map
%input_swtch = 2; % Pick x_section by entering coordinates

load CCP_5nodese.mat
%%%%%%%%%%Map for picking x_section location%%%%%%%%%%%%%%
ilatlim = [(min(lat))-.1 (max(lat))+.1];
ilonlim = [(min(lon))-.1 (max(lon))+.1];
ind = find(z == 20);
if input_swtch == 1;
  base_map_v2(ilatlim, ilonlim, Stns,xlpierce(ind,:),ylpierce(ind,:));
%%%%%%%Pick Locations for x_section ends%%%%%%%%%%
  [yy xx] = ginput;
  close all;
%%%%%Enter coordinate  for cross section here %%%%%%%
%% 
elseif input_swtch == 2;
  xx = [-85.8 -85.06]'; %x coordinates
  yy = [9.98 10.72]'; %y coordinates
end
%%%%%%%Plot X_section Location%%%%%%%%%%%%%%5
  base_map_v2(ilatlim, ilonlim, Stns,xlpierce(ind,:),ylpierce(ind,:))
  geoplot(yy,xx,'Color','blue','LineWidth',5)
  %geobasemap('colorterrain');
%%%%%%%load RF Grid%%%%%%%%%%%
  z = -z;
  [xrf,yrf,zrf]=meshgrid(lon,lat,z);
  cmp = (permute(CCP,[2,3,1]));
  cmp_hits = (permute(HITS,[2,3,1]));
  bstd = (permute(BSTD,[2,3,1]));

%%%%%%%%Create Surface for Cross Section%%%%%%%%%%%%
x_val = []; y_val = []; d = [];
for i = 1:(size(xx,1)-1);
 d = [d distance(yy(i),xx(i),yy(i+1),xx(i+1))];
end
for i = 1:(size(xx,1)-1);
        x_val = [x_val linspace(xx(i),xx(i+1),1000*round(d(i)./sum(d),3))];
        y_val = [y_val linspace(yy(i),yy(i+1),1000*round(d(i)./sum(d),3))];
end
rf_dist = 0;
for i = 2:length(x_val)
        rf_dist(i) = lldistkm([y_val(i-1),x_val(i-1)],[y_val(i),x_val(i)]);
end
rf_dist = cumsum(rf_dist);
td = (min(z):1:0)';

        for i = 1:size(x_val,2);
          zd(:,i) = td;
        end;
        for i = 1:size(td,1);
          xd(i,:) = x_val';
          rxd(i,:) = rf_dist';
          yd(i,:) = y_val';
          ryd(i,:) = rf_dist';
        end;
        clear td

%%%%%%%%Project RFs to Cross Section Surface%%%%%%%%%%%
        rf_angle = -interp3(xrf,yrf,zrf,cmp,xd,yd,zd); %
        rf_hits_angle = interp3(xrf,yrf,zrf,cmp_hits,xd,yd,zd);
        rf_bstd_angle = interp3(xrf,yrf,zrf,bstd,xd,yd,zd);
%%%%%%%Vertically Smooth RF Cross Section%%%%%%%%%%
        yn = find(isnan(rf_angle) == 1);
        for kk = 1:size(rf_angle,2); 
           rf_angle(:,kk) = smooth(rf_angle(:,kk),3);
        end
        rf_angle(yn) = NaN; 
        topo = topography_angle_v4(x_val',y_val',rf_dist,'default','noplot');
        moho = moho2_dep(x_val',y_val','Moho.mat');
        lvz = moho2_dep(x_val',y_val','LVZ.mat');
%%%%%%Plot RFs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
      amp_scl = 40;%scale for RF trace amps
      hold on
         p1 = pcolor(rf_dist,(zd(:,1)),rf_angle); set(p1,'edgecolor','none');
        for k = 1:50:1000
          xda(1:(size(rf_angle,1)),1) = rf_dist(1,k);
          xamp = xda + amp_scl*rf_angle(:,k);
          bmxamp = xda + amp_scl*(rf_angle(:,k)+rf_bstd_angle(:,k));
          bmnamp = xda + amp_scl*(rf_angle(:,k)-rf_bstd_angle(:,k));
          plot(xamp,(zd(:,1)),'k')
          plot(bmxamp,(zd(:,1)),'k:','LineWidth',0.25)
          plot(bmnamp,(zd(:,1)),'k:','LineWidth',0.25)
        end
  %%%%%%%%%Plot Peaks%%%%%%      
       %rf_peaks(rf_angle,rf_dist,zd,20,55); 
%%%%%%%%%%%Set Color Map%%%%%%%%
       load NaNmap.mat
       load No_Green
       cprf = no_green;
%       colormap(NaNmap); sca = 0.30; caxis([-sca sca]);
       colormap(cprf); sca = 0.5; caxis([-sca sca]);  
%%%%%%%%%%Draw in Topo%%%%%%%%%%
        topo = topography_angle_v4(x_val',y_val',rf_dist,'default');
        plot(rf_dist',-moho,'r--','LineWidth',4)
        plot(rf_dist',-lvz,'b--','LineWidth',4)
        [bthick] = basement_multiples(rf_dist,x_val,y_val,min(lon),max(lon),min(lat),max(lat));
        %Stn_Plot_v2(x_val',y_val',rf_dist,Stns,bin_sz*bin_sm)
        %EQ_Plot_Mat(x_val',y_val',rf_dist,bin_sz*bin_sm,'Event_Data.mat')
        %EQ_Plot_v3_fetch('2010-01-01 00:30:00','2020-3-31 23:30:00',5,[ilatlim ilonlim],'ISC',x_val',y_val',rf_dist,bin_sz*bin_sm*3);
%%%%%%%%%Frame and Adjust Figure%%%%%%%%%
        colorbar
        daspect([1,1,1]);
        xlim([min(rf_dist),max(rf_dist)]);
        box on
        ylim([min(z),40])
        daspect([1,1,1])
        xlabel('Distance (km)')
        set(gca,'FontSize',18)
        ylabel('Depth (km)')
        yticks([min(z):10:0])
        daspect([1,1,1]);
%%%%%%Plot RF HITS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
      hold on
       p1 = pcolor(rf_dist,(zd(:,1)),rf_hits_angle); set(p1,'edgecolor','none');
%%%%%%%%%%%Set Color Map%%%%%%%%
       load NaNmap.mat
       load No_Green
       cprf = no_green;
%       colormap(NaNmap); sca = 0.30; caxis([-sca sca]);
       colormap(flipud(hot)); sca = 500; caxis([0 sca]);
%%%%%%%%%%Draw in Topo%%%%%%%%%%
        topo = topography_angle_v4(x_val',y_val',rf_dist,'default');
        Stn_Plot_v2(x_val',y_val',rf_dist,Stns,bin_sz*bin_sm)
%%%%%%%%%Frame and Adjust Figure%%%%%%%%%
        colorbar
        daspect([2,1,1]);
        xlim([min(rf_dist),max(rf_dist)]);
        box on
        %ylim([min(z),50])
        ylim([min(z),30])
        

%%%%%%Plot RF Bootstrap%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
      hold on
      p1 = pcolor(rf_dist,(zd(:,1)),rf_bstd_angle); set(p1,'edgecolor','none');

%%%%%%%%%%%Set Color Map%%%%%%%%
       load NaNmap.mat
       load No_Green
       cprf = no_green;
%       colormap(NaNmap); sca = 0.30; caxis([-sca sca]);
       colormap(flipud(hot)); sca = .2; caxis([0 sca]);
%%%%%%%%%%Draw in Topo%%%%%%%%%%
        topo = topography_angle_v4(x_val',y_val',rf_dist,'default');
        Stn_Plot_v2(x_val',y_val',rf_dist,Stns,bin_sz*bin_sm)
        %EQ_Plot_Mat(xd,yd,x_min,x_max,y_min,y_max,'EQs.txt',bin_sz*bin_sm)
%%%%%%%%%Frame and Adjust Figure%%%%%%%%%
        colorbar
        daspect([2,1,1]);
        xlim([min(rf_dist),max(rf_dist)]);
        box on
        ylim([min(z),30])


