function [] = base_map(ilatlim, ilonlim, Stns,xpierce,ypierce)
  gx = geoaxes;
  geobasemap(gx,'topographic')
  geolimits(gx,ilatlim,ilonlim)
  set(gx,'Grid','on');
  hold on
%%%%%%%Plot topgraphy cities and volcanoes%%%%%%%%
  %C = surfm(ilatlim, ilonlim,SA);
  cities = shaperead('worldcities', 'UseGeoCoords', true);
  geoscatter(vertcat(cities.Lat),vertcat(cities.Lon),'Marker', 'o', 'MarkerFaceColor', 'blue','MarkerEdgeColor','k')
  load volcano.txt
  geoscatter(volcano(:,1), volcano(:,2),'Marker','^','MarkerFaceColor', 'red','MarkerEdgeColor','k')
%%%%%%%Plot Pierce Points%%%%%%%%
if ~isempty(xpierce)
    geoscatter(ypierce, xpierce,'Marker', 'x','MarkerFaceColor', 'r','MarkerEdgeColor','k')
end
%%%%%%%Plot Stations%%%%%%%%
 stn_names = (fields(Stns));
   for i = 1:length(stn_names);
     stn_name = char(stn_names(i));
     stn_lat(i) = Stns.(stn_name).Station_Data.Latitude;
     stn_lon(i) = Stns.(stn_name).Station_Data.Longitude;
   end
   if exist('stn_lat','var')
   geoscatter(stn_lat, stn_lon,'Marker', 'd','MarkerFaceColor', 'g','MarkerEdgeColor','k')
   end
%%%%%%%%%END%%%%% 

%%% Plot Whitmyere Lines%%%
%S = shaperead('./Data/Whitmeyer2007_Units_20080822.dbf');

%for i = 1:size(S,1)
%    xx = S((i)).X;
%    yy = S((i)).Y;
%    zone = repmat(15,size(xx));
%    [lat lon] = whit_utm2deg(xx,yy,(zone));
%    geoplot(lat,lon,'k') 
%end
%%%%%end%%%

