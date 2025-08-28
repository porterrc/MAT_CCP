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
         plot(evt_list(:,1),(evt_list(:,2)),'k^','MarkerFaceColor','c','MarkerSize',6);
         ht = text(evt_list(:,1),(evt_list(:,2)+.5),names,'Interpreter', 'none');
    end
    set(ht,'Rotation',45)
%%%%%%%END Plot EVENTS%%%%%%%%%%%
