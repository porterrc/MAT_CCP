function [] = plot_array(X,Y,z,stnx,stny,sname,xpierce,ypierce,sx,sy,sz,svmodel,pdepth)
figure
hold on

if size(svmodel,1) == length(sx);
    [SY SX] = meshgrid(sy,sx);
elseif size(svmodel,1) == length(sy);
    [SX SY] = meshgrid(sx,sy);
end
ind = min(find(abs(sz-pdepth) == min(abs(sz-pdepth))));
c = pcolor(SX,SY,squeeze(svmodel(:,:,ind)));
set(c,'Edgecolor','none');
plot(X,Y,'+k');
plot(stnx,stny,'*');
text(stnx,stny,sname);
ind = find(z == pdepth);
plot(xpierce(ind,:),ypierce(ind,:),'r.');
colorbar
set(gca,'Zdir','reverse')