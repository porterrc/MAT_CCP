close all; clear all;
addpath ~/Dropbox/Data/Topo_Files
addpath ./Functions
addpath ./Data
addpath ./Functions/Color_Palette

decon = 'ID' % Choose Deconvolution method

lat1 = 34.5; % southern lat for grid
lat2 = 36.4 ;% northern lat for grid
lon1 = -112.7; % western edge for grid 
lon2 = -110.6 ;% eastern edge for grid
input_swtch = 2;
sf = 8; %rf amp scaling factor
min_amp = -2;
max_amp = 2;
load RF_5_Merge.mat


%%%%%%%%%%Map for picking x_section location%%%%%%%%%%%%%%
ilatlim = [(min(lat1))-.1 (max(lat2))+.1];
ilonlim = [(min(lon1))-.1 (max(lon2))+.1];
if input_swtch == 1;
  base_map_v2(ilatlim, ilonlim, Stns,[],[]);
%%%%%%%Pick Locations for x_sectcd ~/ion ends%%%%%%%%%%
  [yy xx] = ginput;
  close all;
%%%%%Enter coordinate  for cross section here %%%%%%%
elseif input_swtch == 2;
  %xx = [ -112.33 -111.19]'; %x coordinates
  %yy = [35.17 35.435]'; %y coordinates   
   xx = [ -112.521 -110.916]'; %x coordinates
   yy = [35.132 35.5]'; %y coordinates
end
%%%%%%%Plot X_section Location%%%%%%%%%%%%%%5
  base_map_v2(ilatlim, ilonlim, Stns,[],[])
  geoplot(yy,xx,'Color','blue','LineWidth',5)
  %geobasemap('colorterrain');

%%%%%%%%Create Surface for Cross Section%%%%%%%%%%%%
x_val = []; y_val = []; d = [];
for i = 1:(size(xx)-1);
 d = [d distance(yy(i),xx(i),yy(i+1),xx(i+1))];
end
for i = 1:(size(xx)-1);
        x_val = [x_val linspace(xx(i),xx(i+1),1000*round(d(i)./sum(d),3))];
        y_val = [y_val linspace(yy(i),yy(i+1),1000*round(d(i)./sum(d),3))];
end
rf_dist = 0;
for i = 2:length(x_val)
        rf_dist(i) = lldistkm([y_val(i-1),x_val(i-1)],[y_val(i),x_val(i)]);
end
rf_dist = cumsum(rf_dist);

[S_Line S_dist] = Stn_Line(x_val,y_val,rf_dist,Stns,10);

figure
hold on
for i = 1:length(S_Line);
    rf = -trimmean(Data.(char(S_Line(i))).RFs_in_Time.(decon).RFs,10);
    tvec = Data.(char(S_Line(i))).RFs_in_Time.(decon).Time';
    ind = find(tvec >= -0.2);
    rf = rf(ind); tvec = tvec(ind);
    if max(rf) > max_amp | min(rf) < min_amp
        continue
    end
    % Extract positive and negative part
    rfp = (rf + abs(rf))/2;
    rfp(1) = 0; rfp(end) = 0;
    rfn = (rf - abs(rf))/2;
    rfn(1) = 0; rfn(end) = 0;
    patch(S_dist(i)+rfp*sf,tvec,'red')
    patch(S_dist(i)+rfn*sf,tvec,'blue')
end
 set(gca,'Ydir','reverse')
 basement_multiples_time(Stns,S_Line,S_dist,4,2.29,'k');
 basement_multiples_time(Stns,S_Line,S_dist,3,1.71,[.7 .7 .7]);
topo = topography_angle_v4(x_val',y_val',rf_dist,'default');
Stn_Plot(x_val',y_val',rf_dist,Stns,10)
ylim([-5 10])
xlim([0 max(rf_dist)])

box on
set(gca,'FontSize',32)
xlabel('Distance (km)')
ylabel('Time (s)')
yticks([0:2:10])


function [topokm] = topography_angle_v4(x_val,y_val,rf_dist,resolution,dplot);
%[SA,refvec] = gtopo30('~/Dropbox/Data/Topo_Files/',2,[min(y_val)-.5,max(y_val)+.5],[min(x_val)-.5,max(x_val)+.5]);
[SA,lon,lat] = TopoDownload([min(x_val)-.5,max(x_val)+.5],[min(y_val)-.5,max(y_val)+.5],resolution);

i = find(isnan(SA));
SA(i)=zeros(size(i));

sz_SA = size(SA);
X = lon;
Y = lat;

topo = -interp2(X',Y',SA',x_val,y_val);
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
end

function [] = Stn_Plot(x_val,y_val,rf_dist,Stns,binsz)
 stn_names = (fields(Stns));
 j = 1;
 for i = 1:length(stn_names);
   stn_name = [char(stn_names(i))];
   name(i,:) = [cell(stn_names(i))];
   stn_lat(i) = Stns.(stn_name).Station_Data.Latitude;
   stn_lon(i) = Stns.(stn_name).Station_Data.Longitude;
   stn_elv(i) = Stns.(stn_name).Station_Data.Elevation;
end
   eq = [stn_lon' stn_lat' stn_elv']; 
    evt_list = [];
    x = x_val; y = y_val;
    [idx d] = knnsearch([x y],[eq(:,1) eq(:,2)]);
    i = find(d < binsz/111);
    idx = idx(i); d = d(i); z = eq(i,3); names=name(i); 
    evt_list = [x(idx) y(idx) z/1000];
    evt_list = [rf_dist(idx)' z/1000];
    if size(evt_list,1) >= 1;
         plot(evt_list(:,1),-(evt_list(:,2)),'k^','MarkerFaceColor','b','MarkerSize',12);
         %ht = text(evt_list(:,1),(evt_list(:,2)+.5),names,'Interpreter', 'none');
    end
    %set(ht,'Rotation',45)
%%%%%%%END Plot EVENTS%%%%%%%%%%%
end