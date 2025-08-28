function [slabd] = slab2_dep(x_val,y_val)
    slab = load('/Users/rp522/Dropbox/Data/Topo_Data/Topo_Files/sam_slab2_dep_02.23.18.xyz');
    x = slab(:,1)-360;
    y = slab(:,2);
    z = slab(:,3);
    xd = length(unique(x));
    yd = length(unique(y));
    X = reshape(x,[xd yd])'; 
    Y = reshape(y,[xd yd])';
    Z = reshape(z,[xd yd])';
    
   slabd = interp2(X,Y,Z,x_val,y_val);