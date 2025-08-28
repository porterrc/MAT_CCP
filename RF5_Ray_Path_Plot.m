close all; clear all;
addpath ./Functions
addpath ./Functions/irisFetch-matlab-2.0.12
addpath ./Data
addpath ./Functions/Color_Palette
%
%input_swtch = 1; % Pick x_section location on map
input_swtch = 2; % Pick x_section by entering coordinates

load CCP_ray_info_5s.mat
binsz = 20

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
  xx = [ -112.33 -111.19]'; %x coordinates
  yy = [35.17 35.435]'; %y coordinates
end
%%%%%%%Plot X_section Location%%%%%%%%%%%%%%
  figure
  base_map_v2(ilatlim, ilonlim, Stns,xlpierce(ind,:),ylpierce(ind,:))
  geoplot(yy,xx,'Color','blue','LineWidth',5)
  %geobasemap('colorterrain');
%%%%%%%load RF Grid%%%%%%%%%%%
  z = -z;
  [xrf,yrf,zrf]=meshgrid(lon,lat,z);

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

        figure
        hold on
        Ray_Path_Plot_Amp(x_val,y_val,rf_dist,binsz,xlpierce,ylpierce,z,seis);         
        topo = topography_angle_v4(x_val',y_val',rf_dist,'default');
        %Stn_Plot_v2(x_val',y_val',rf_dist,Stns,10)
        xlim([0 max(rf_dist)])
        ylim([-80 50])


function [] = Ray_Path_Plot_Amp(x_val,y_val,rf_dist,binsz,xpierce,ypierce,z,seis)
    rp_plot = nan(size(xpierce));   
    for i = 1:length(z)
        rp = [xpierce(i,:)'  ypierce(i,:)'];
        samp = seis(i,:);
        evt_list = [];
        x = x_val'; y = y_val';
        [idx d] = knnsearch([x y],[rp(:,1) rp(:,2)]);
        ind = find(d >= binsz/111);
        ind0  = find(isnan(d));
        %idx(ind) = nan;
        rp_plot(i,:) = rf_dist(idx);
        rp_seis(i,:) = samp(idx);
        rp_plot(i,ind) = nan;
        rp_plot(i,ind0) = nan;
        
    end
         load No_Green
         cprf = no_green;
         %scatter(rp_plot,z'.*ones(size(rp_plot)),10,rp_seis,'filled')
         s = scatter(reshape(rp_plot,[],1),reshape(z'.*ones(size(rp_plot)),[],1),30,(reshape(-(rp_seis),[],1)),"filled")
         s.MarkerFaceAlpha = .4;
         colormap(cprf); sca = 1; caxis([-sca sca]);  
    end
%%%%%%%END Plot EVENTS%%%%%%%%%%%