[CCPiso CCPan CCPaz] = CCP_anistroy(CCP,CCP_RAW,BAZ_RAW)
CCPan = nan(size(CCP)); CCPiso = nan(size(CCP)); CCPaz = nan(size(CCP));
xrng = size(CCP_RAW,2);
yrng = size(CCP_RAW,3);
parfor zz=1:size(CCP_RAW,1)
    for xx = 1:xrng
        for yy = 1:yrng
            cellamp = squeeze(CCP_RAW(zz,xx,yy,:));
            cellbaz = squeeze(BAZ_RAW(zz,xx,yy,:));
            if all(isnan(cellamp)) || length(cellamp(~isnan(cellamp))) < 3
                continue
            end
            sfit = sinFit(cellbaz,cellamp);
            CCPan(zz,xx,yy) = sfit.a;
            CCPaz(zz,xx,yy) = sfit.b;
            CCPiso(zz,xx,yy) = mean(cellamp-sfit.a.*(sind(cellbaz+sfit.b*180/pi)),'omitnan');
        end
    end
end