function [pks,ws,ps,dpk] = Moho_Picker(CCP,lon,lat,z,mnd,mxd);
zd = z;
pks = nan(size(CCP,[2:3]));
ws = nan(size(CCP,[2:3]));
pr = nan(size(CCP,[2:3]));
for i = 1:size(CCP,2);
    ilat = lat(i);
    for j = 1:size(CCP,3)
        ilon = lon(j);
        d = find(zd >= mnd & zd <= mxd);
        rfd = CCP(:,i,j);
        nanind = isnan(rfd);
        rfd = smooth(rfd,9);
        rfd(nanind) = nan;
        rfd = rfd(d);
        if min(isnan(rfd)) == 0
            [ps,ls,w,p] = findpeaks(rfd,'MinPeakHeight',.075,'MinPeakDistance',5);
            if(size(p,1) > 0)
                ind = find(ps == max(ps));
                dpk(i,j) = (z(d(ls(ind))));
                pks(d(ls),i,j) = ps;
                ws(d(ls),i,j) = w;
                pr(d(ls),i,j) = p;
            else
                dpk(i,j) = nan;
                pks(:,i,j) = nan;
            end
        else
            dpk(i,j) = nan;
            pks(:,i,j) = nan;
        end
    end
end
