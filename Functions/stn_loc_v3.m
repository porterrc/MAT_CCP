function [stn_named sx sy se] = stn_loc(Stns,lat1,lat2,lon1,lon2,dd)
 stn_names = (fields(Stns));
 j = 1;
 for i = 1:length(stn_names);
   stn_name = char(stn_names(i));
   stn_lat = Stns.(stn_name).Station_Data.Latitude;
   stn_lon = Stns.(stn_name).Station_Data.Longitude;
   stn_elv = Stns.(stn_name).Station_Data.Elevation;
   if stn_lon < lon2 & stn_lon > lon1 & stn_lat > lat1 & stn_lat < lat2; 
     sy(j) = deg2km(stn_lat-lat1);
     sx(j) = dd*(stn_lon-lon1);
     se(j) = stn_elv/1000;
     stn_named(j) = stn_names(i);
     j = j + 1;
   end
end
