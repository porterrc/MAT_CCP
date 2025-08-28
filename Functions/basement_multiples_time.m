function [thick] = basement_multiples_time(Stns,S_Line,S_Dist,Vpb,Vsb,lc);


 Stk = ncread('Pcb.nc','PcB');
 Slon = ncread('Pcb.nc','lon');
 Slat= ncread('Pcb.nc','lat'); 
 Stk(Stk == -9999) = nan;
 [Sy Sx] = meshgrid(Slat,Slon);
  

for i = 1:length(S_Line);
    selv(i) = Stns.(char(S_Line(i))).Station_Data.Elevation;
    xq(i) = Stns.(char(S_Line(i))).Station_Data.Longitude;
    yq(i) = Stns.(char(S_Line(i))).Station_Data.Latitude;
end

   St = interp2(Sy,Sx,Stk,yq,xq); 
   thick = (selv-St*0.3048)/1000;

 k = find(thick < 0); thick(k) = 0; % remove sed thicknesses less than 0
 p = 0.06; % assume vertical incidence
 t1a = thick*((1/Vsb^2 -p^2)^.5-(1/Vpb^2 -p^2)^.5); %time of ps
 t1b = thick.*((1/Vsb^2 -p^2)^.5+(1/Vpb^2 -p^2)^.5); % timing of 1st multiple
 t2b = 2.*thick.*(1/Vsb^2 -p^2)^.5; % timing of 2nd  multiple

 [~,s_ind] = sort(S_Dist);
plot(S_Dist(s_ind),t1a(s_ind)','_','Color',lc,'LineWidth',4)
plot(S_Dist(s_ind),t1b(s_ind)',':','Color',lc,'LineWidth',4)
plot(S_Dist(s_ind),t2b(s_ind)','--','Color',lc,'LineWidth',4)

