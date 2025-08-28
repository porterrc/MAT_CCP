function [sx sy sz svmodel] = shear_loc(file,dd)
     load file
     ilat = ncread(file,'latitude');
     ilon = ncread(file,'longitude');
     svmodel = ncread(file,'vs');
     sy = deg2km(ilat-lat1);
     sx = dd*(ilon-lon1);
     sz = ncread(file,'depth');
