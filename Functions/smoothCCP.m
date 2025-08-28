function [CCP] = smoothCCP(CCP,minwidth,maxwidth,dzi);

for i = 1:size(CCP,2)
    i
    for j = 1:size(CCP,3)
        dvec = squeeze((CCP(:,i,j)));
        ind = isnan(dvec);
        dvec(ind) = 0;
        dvec = bandpass(dvec,[1/maxwidth,1/minwidth],dzi);
        dev(ind) = nan;
        CCP(:,i,j) = dvec; 
    end
end