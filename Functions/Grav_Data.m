function [X,Y,Z] = Grav_Data(file)

xrng = ncread(file,'x_range');
yrng = ncread(file,'y_range');
inc = ncread(file,'spacing');
z = ncread(file,'z');

x = xrng(1):inc(1):xrng(2);
y = yrng(2):-inc(2):yrng(1);

[Y,X] = meshgrid(y,x);
Z = reshape(z,size(X));