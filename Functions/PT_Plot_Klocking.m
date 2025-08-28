function [] = PT_Plot(mmp,x,y,z,xd,yd,zd,rf_dist,file,binsz)
    pt = load(file);
    xpt = xd(1,:)'; ypt = yd(1,:)'; % define x and y pts
    
    
    evt_list = [];
    [idx d] = knnsearch([xpt ypt],[pt(:,2) pt(:,1)]); % find distances between samples and transect
    i = find(d < binsz/111);
    idx = idx(i); %choose samples within range
    d = d(i); % choose only distances within range
    
    zp = -mean([pt(i,3),pt(i,4)],2);
    t = mean([pt(i,5),pt(i,6)],2)+0.48.*-zp; %temp
    
    
        %find PT for error bars    % fine min and max temps for error bars
    mnt = pt(i,5)+0.48.*-zp; 
    mxt = pt(i,6)+0.48.*-zp; % min and max temp of calcs
    mnz = -pt(i,3); % mn depth
    mxz = -pt(i,4); % mn depth
    
    %%%Calc Error Bar Lengths
    yneg = zp-mnz; ypos = mxz-zp;
    xneg = t-mnt; xpos = mxt-t;
    xpos = xpos/2000*(max(rf_dist)); 
    xneg = xneg/2000*(max(rf_dist));
    
    %make matrix of results
    evt_list = [rf_dist(idx)' zp t]; 
        if size(evt_list,1) >= 1;
           %errorbar(evt_list(:,1),evt_list(:,2),yneg,ypos,xneg,xpos,'k.','CapSize',10)
           errorbar(evt_list(:,1),evt_list(:,2),yneg,ypos,'k.','CapSize',10)
           hold on
           scatter(evt_list(:,1),(evt_list(:,2)),60,evt_list(:,3),'s','filled','markeredgecolor','black');
           %errorbar(max(rf_dist)/10,(min(ylim)+ 20),50/2000*(max(rf_dist)),'horizontal','k')
           %text(max(rf_dist)/10,(min(ylim)+ 14),['50' char(176) ' C'], 'HorizontalAlignment','Center','FontSize',16);
        end

