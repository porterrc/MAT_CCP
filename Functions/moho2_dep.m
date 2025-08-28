function [slabd slaba] = moho2_dep(x_val,y_val,file)
    load(file);
    x = Lon;
    y = Lat;
    z = idepth;
    za = iamp;
    
   slabd = interp2(x,y,z,x_val,y_val);
   slaba = interp2(x,y,za,x_val,y_val);