function [sx sy sz svmodel] = shear_loc(file,lat1,lat2,lon1,lon2,dd)
     ilat = ncread(file,'latitude');
     ilon = ncread(file,'longitude');
     svmodel = ncread(file,'vs');
     svmodel(find(svmodel == 9999)) = nan;
     sy = deg2km(ilat-lat1);
     sx = dd*(ilon-lon1);
     sz = ncread(file,'depth');
