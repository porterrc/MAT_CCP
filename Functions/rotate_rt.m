function [r,t] = rotate_rt(n,e,phi,tsnaz) 
    phi = phi - tsnaz;
    cphi = cosd(phi);
    sphi = sind(phi);
    
    r = -(cphi*n+sphi*e);
    t = -sphi*n+cphi*e;
end