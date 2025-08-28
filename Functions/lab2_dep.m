function [slabd] = lab2_dep(x_val,y_val)
    slab = load('./Data/Thickness210e3_75.txt');
    x = slab(:,1);
    y = slab(:,2);
    z = slab(:,3);
    xd = length(unique(x));
    yd = length(unique(y));
    X = reshape(x,[yd xd])'; 
    Y = reshape(y,[yd xd])';
    Z = reshape(z,[yd xd])';
    
   slabd = interp2(X',Y',Z',x_val,y_val);