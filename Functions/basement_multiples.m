function [thick] = basement_multiples(rf_dist,xd,yd,x_min,x_max,y_min,y_max);

% Define topo
[SA,lon,lat] = TopoDownload([min(xd)-.5,max(xd)+.5],[min(yd)-.5,max(yd)+.5],'default');
i = find(isnan(SA));
SA(i)=zeros(size(i));

sz_SA = size(SA);
X = lon;
Y = lat;
topo = interp2(X',Y',SA',xd,yd);
a = size(topo);
slength = round(60/(rf_dist(2)-rf_dist(1))); % smooths to 60 km
topokm = smooth(topo,slength);


 Stk = ncread('Pcb.nc','PcB');
 Slon = ncread('Pcb.nc','lon');
 Slat= ncread('Pcb.nc','lat'); 
 Stk(Stk == -9999) = nan;

 xq = (x_min-.5):.01:(x_max+.5);
 yq = (y_min-.5):.01:(y_max+.5);
 [xq,yq] = meshgrid(xq,yq);


[Sy Sx] = meshgrid(Slat,Slon);

 St = interp2(Sy,Sx,Stk,yq,xq);

 Sl = [reshape(xq,[numel(xq),1]) reshape(yq,numel(yq),1)];
 St = reshape(St,[numel(St),1]); 
 xl = [xd(1,:)' yd(1,:)'];
 k = dsearchn(Sl,xl);
 elev = St(k); %thickness value for evey x y point
 thick = (topokm'-elev'*0.3048)/1000;
 k = find(thick < 0); thick(k) = 0; % remove sed thicknesses less than 0
 Vpb = 3.5; % basin Vp
 Vsb = 1.9; % basin Vs
 p = 0; % assume vertical incidence
 t1a = thick*((1/Vsb^2 -p^2)^.5-(1/Vpb^2 -p^2)^.5); %time of ps
 t1b = thick.*((1/Vsb^2 -p^2)^.5+(1/Vpb^2 -p^2)^.5); % timing of 1st multiple
 t2b = 2.*thick.*(1/Vsb^2 -p^2)^.5; % timing of 2nd  multiple

% read in vel model to migrate time to depth
 vmodel = 'iasp91.tvel';
 dzi = 0.5;
 z_max = 100;
 vel = dlmread(vmodel,'', 2, 0); % read in velocity mode
 z = 0:dzi:z_max; % create depth vector
 tmp = vel(:,1);
 [~,ind,~] = unique(tmp); tmp(ind) = -1; ind = find(tmp > 0); % account for duplicate values
 vel(ind,1) = vel(ind,1)+0.0001; % account for duplicate values
 vp = interp1(vel(:,1),vel(:,2),z(1:end)+0.5*dzi); % interpolate vp and vs to depth
 vs = interp1(vel(:,1),vel(:,3),z(1:end)+0.5*dzi); % interpolate vp and vs to depth
 qa = sqrt(1./vp.^2 - p.^2);
 qb = sqrt(1./vs.^2 - p.^2);
 dt = dzi*(qb - qa); % time difference at each depth interval
 tbot = cumsum(dt); % get the time at bottom of each layer
 ttop = tbot - dt;  %get time at top of each layer

for i = 1:size(t1b,2);
 k = find(t1a(i)<tbot);
 d0(i) = dzi*k(1);
 k = find(t1b(i)<tbot);
 d1(i) = dzi*k(1);
 k = find(t2b(i)<tbot);
 d2(i) = dzi*k(1);
end

d0 = d0-topokm'./1000;
d1 = d1-topokm'./1000;
d2 = d2-topokm'./1000;

plot(rf_dist,-smooth(d0,11),'r_','LineWidth',2)
plot(rf_dist,-smooth(d1,11),'r:','LineWidth',2)
plot(rf_dist,-smooth(d2,11),'b:','LineWidth',2)

