
close all; clear all;
addpath ./Functions
addpath ./Functions/irisFetch-matlab-2.0.12
addpath ./Data
addpath ./Functions/Color_Palette
%
%input_swtch = 1; % Pick x_section location on map
input_swtch = 2; % Pick x_section by entering coordinates


load CCP_5nodes.mat
slab = load('/Users/rp522/Dropbox/Data/Topo_Data/Topo_Files/sam_slab2_dep_02.23.18.xyz');
shear_file = 'WUS324.r0.0.nc'

%%%%%%%%%%Map for picking x_section location%%%%%%%%%%%%%%
ilatlim = [(min(lat))-.1 (max(lat))+.1];
ilonlim = [(min(lon))-.1 (max(lon))+.1];
ind = find(z == 20);
if input_swtch == 1;
  figure
  base_map_v2(ilatlim, ilonlim, Stns,xlpierce(ind,:),ylpierce(ind,:));
%%%%%%%Pick Locations for x_section ends%%%%%%%%%%
  [yy xx] = ginput;
  close all;
%%%%%Enter coordinate  for cross section here %%%%%%%
elseif input_swtch == 2;
  xx = [ -112.521 -110.916]'; %x coordinates
  yy = [35.132 35.5]'; %y coordinates
end
%%%%%%%Plot X_section Location%%%%%%%%%%%%%%
  figure
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
        x_val = [x_val linspace(xx(i),xx(i+1),2500*round(d(i)./sum(d),3))];
        y_val = [y_val linspace(yy(i),yy(i+1),2500*round(d(i)./sum(d),3))];
end
rf_dist = 0;
for i = 2:length(x_val)
        rf_dist(i) = lldistkm([y_val(i-1),x_val(i-1)],[y_val(i),x_val(i)]);
end
rf_dist = cumsum(rf_dist);
td = (min(z):.25:0)';

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
        rf_angle = interp3(xrf,yrf,zrf,cmp,xd,yd,zd);
        rf_hits_angle = interp3(xrf,yrf,zrf,cmp_hits,xd,yd,zd);
        rf_bstd_angle = interp3(xrf,yrf,zrf,bstd,xd,yd,zd);
%%%%%%%Vertically Smooth RF Cross Section%%%%%%%%%%
        yn = find(isnan(rf_angle) == 1);
        for kk = 1:size(rf_angle,2); 
           rf_angle(:,kk) = -smooth(rf_angle(:,kk),3);
           %rf_angle(:,kk) = -(rf_angle(:,kk));
        end
        rf_angle(yn) = NaN; 
        topo = topography_angle_v4(x_val',y_val',rf_dist,'default','noplot');
        slab = slab2_dep(x_val',y_val');

%%%%%Read Shear File
     svlat = ncread(shear_file,'latitude');
     svlon = ncread(shear_file,'longitude');
     svmodel = permute(ncread(shear_file,'VS'),[1 2 3]);
     pvmodel = permute(ncread(shear_file,'VP'),[1 2 3]);
     svz = double(-ncread(shear_file,'depth'));
     [xsv,ysv,zsv] = meshgrid(svlon,svlat,svz);
     if max(size(xsv) ~= size(svmodel))
         svmodel = permute(svmodel,[2,1,3]);
         pvmodel = permute(pvmodel,[2,1,3]);         
     end

%%%%%%%%Project SV to Cross Section Surface and Smooth%%%%%%%%%%%
    sv_angle = interp3(xsv,ysv,zsv,svmodel,xd,yd,zd);
    pv_angle = interp3(xsv,ysv,zsv,pvmodel,xd,yd,zd);
%%%%%%Plot RFs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
      amp_scl = 20;%scale for RF trace amps
      hold on
         p1 = pcolor(rf_dist,(zd(:,1)),rf_angle); set(p1,'edgecolor','none');
        for k = 1:250:2500
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
       colormap(cprf); sca = 0.4; caxis([-sca sca]);  
%%%%%%%%%%Draw in Topo%%%%%%%%%%
        topo = topography_angle_v4(x_val',y_val',rf_dist,'default');
        plot(rf_dist',slab,'k','LineWidth',4);
        [bthick] = basement_multiples_shear(rf_dist,x_val,y_val,zd,min(lon),max(lon),min(lat),max(lat),sv_angle,pv_angle);
        %plot(rf_dist',slab,'k','LineWidth',4)
        %Stn_Plot_v2(x_val',y_val',rf_dist,Stns,bin_sz*bin_sm)
        %EQ_Plot_Mat(x_val',y_val',rf_dist,bin_sz*bin_sm,'Event_Data.mat')
        EQ_Plot_Mat(x_val,y_val,rf_dist,10,'Sunset_Swarm.mat')
        EQ_Plot_Text(x_val,y_val,rf_dist,10,'SC_EQ_cat.txt')
        PT_Plot_v2(45,x,y,z,xd,yd,zd,rf_dist,'Melt_PT.mat',25,1.5)
        PT_Plot_Plank(45,x,y,z,xd,yd,zd,rf_dist,'Plank_PT.txt',25)
        PT_Plot_Lexi(45,x,y,z,xd,yd,zd,rf_dist,'Lexi_PT.txt',25)
        PT_Plot_Klocking(45,x,y,z,xd,yd,zd,rf_dist,'Klocking_PT.txt',25);
        %EQ_Plot_v3_fetch('2010-01-01 00:30:00','2020-3-31 23:30:00',5,[ilatlim ilonlim],'ISC',x_val',y_val',rf_dist,bin_sz*bin_sm*3);
%%%%%%%%%Frame and Adjust Figure%%%%%%%%%
        a = colorbar;
        xlim([min(rf_dist),max(rf_dist)]);        
        %ylim([min(z),40])
        ylim([-80,40])
        box on
        set(gca,'FontSize',18)
        xlabel('Distance (km)')
        ylabel('Depth (km)')
        a.Label.String = 'RF Amplitude';
        yticks([min(z):10:0])
        daspect([1,1,1]);
%%%%%%Plot Shear Velocity%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%Get Color Palette
    load haxby.mat
    tej = (haxby);
 %%%%Plot Cross section   
    figure
    hold on
    p2 = pcolor(rf_dist,(zd(:,1)),sv_angle); set(p2,'edgecolor','none');
    colormap(tej);
    topo = topography_angle_v4(x_val',y_val',rf_dist,'default');
    
    %%% PLot slab and stations
    plot(rf_dist',slab,'k','LineWidth',4)
        %plot(rf_dist',slab,'k','LineWidth',4)
    Stn_Plot_v2(x_val',y_val',rf_dist,Stns,bin_sz*bin_sm)
 %%%%%%%%%Frame and Adjust Figure%%%%%%%%%
    a = colorbar;
    xlim([min(rf_dist),max(rf_dist)]);        
    %ylim([min(z),40])
    ylim([-80,40])
    caxis([min(min(sv_angle)),max(max(sv_angle))]);
    box on
    set(gca,'FontSize',18)
    xlabel('Distance (km)')
    ylabel('Depth (km)')
    a.Label.String = 'S-Velocity (km/s)';
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
%%%%%%%%%%Draw inda Topo%%%%%%%%%%
        topo = topography_angle_v4(x_val',y_val',rf_dist,'default');
        Stn_Plot_v2(x_val',y_val',rf_dist,Stns,bin_sz*bin_sm)
        Ray_Path_Plot(x_val,y_val,rf_dist,10,xlpierce,ylpierce,z)
%%%%%%%%%Frame and Adjust Figure%%%%%%%%%
        a = colorbar;
        xlim([min(rf_dist),max(rf_dist)]);        
        %ylim([min(z),40])
        ylim([-80,40])
        set(gca,'colorscale','log')
        caxis([10 10^5])
        box on
        set(gca,'FontSize',18)
        xlabel('Distance (km)')
        ylabel('Depth (km)')
        a.Label.String = 'Traces per Gridpoint';
        yticks([min(z):10:0])
        daspect([1,1,1]);
        

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
        a = colorbar;
        xlim([min(rf_dist),max(rf_dist)]);        
        ylim([min(z),40])
        %ylim([-80,40])
        box on
        set(gca,'FontSize',18)
        xlabel('Distance (km)')
        ylabel('Depth (km)')
        a.Label.String = 'Bootstrap Variation';
        yticks([min(z):10:0])
        daspect([1,1,1]);

%%%%%%Plot Grid Size%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('GSIZE')
    gsize = (permute(GSIZE,[2,3,1]));
    rf_gz_angle = interp3(xrf,yrf,zrf,gsize,xd,yd,zd);
    figure
      hold on
      p1 = pcolor(rf_dist,(zd(:,1)),rf_gz_angle); set(p1,'edgecolor','none');

%%%%%%%%%%%Set Color Map%%%%%%%%
       load NaNmap.mat
       load No_Green
       cprf = no_green;
%       colormap(NaNmap); sca = 0.30; caxis([-sca sca]);
       colormap(flipud(hot)); sca = 15; caxis([0 sca]);
%%%%%%%%%%Draw in Topo%%%%%%%%%%
        topo = topography_angle_v4(x_val',y_val',rf_dist,'default');
        Stn_Plot_v2(x_val',y_val',rf_dist,Stns,bin_sz*bin_sm)
        %EQ_Plot_Mat(xd,yd,x_min,x_max,y_min,y_max,'EQs.txt',bin_sz*bin_sm)
%%%%%%%%%Frame and Adjust Figure%%%%%%%%%
        a = colorbar;
        xlim([min(rf_dist),max(rf_dist)]);        
        ylim([min(z),40])
        %ylim([-80,40])
        box on
        set(gca,'FontSize',18)
        xlabel('Distance (km)')
        ylabel('Depth (km)')
        a.Label.String = 'Grid Size';
        yticks([min(z):10:0])
        daspect([1,1,1]);
end