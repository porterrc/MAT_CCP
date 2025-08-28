function [] = PT_Plot(ctk,x,y,z,xd,yd,zd,rf_dist,file,binsz,h2o)
    load(file);
    xpt = xd(1,:)'; ypt = yd(1,:)'; % define x and y pts
    
    % choose which water % to plot PT at
    [d1 ~] = find(DDA(:,:,6) == h2o); 
    d1 = unique(d1);
    pt = squeeze(DDA(d1,:,:));
  
    % fine min and max temps for error bars
    mnt = min(DDA(:,:,3)); mxt = max(DDA(:,:,3)); % min and max temp of calcs
    mnp = min(DDA(:,:,4)); mxp = max(DDA(:,:,4)); % min and max temp of calcs
    
    evt_list = [];
    [idx d] = knnsearch([xpt ypt],[pt(:,1) pt(:,2)]); % find distances between samples and transect
    i = find(d < binsz/111);
    idx = idx(i); %choose samples within range
    d = d(i); % choose only distances within range
    pres = nan(size(z));
    pres(z >= -ctk) = (2700*-z(z >= -ctk)*1000*9.81);
    pres(z < -ctk) = (3300*((-z(z < -ctk))-ctk)*1000*9.81)+ctk*2700*1000*9.81;
    zp = zd(:,idx);
    pres = pres';
    
    %calculate pressure in places where not previously constrained
    zp = (-zp)-min(-zp); 
    pres(isnan(pres)) =  zp(isnan(pres)).*1000*2700*9.81;% Pressure in crust
    
    %find closest pressure to PT estimates
    pi = find(abs(pres-pt(i,4)'*10^9) == min(abs(pres-pt(i,4)'*10^9)));
    zp = zd(pi); % depth
    t = pt(i,3); %temp
    
    %find PT for error bars
    pi = find(abs(pres-mnp(i)*10^9) == min(abs(pres-mnp(i)*10^9)));
    mnz = zd(pi); % mn depth
    pi = find(abs(pres-mxp(i)*10^9) == min(abs(pres-mxp(i)*10^9)));
    mxz = zd(pi); % mn depth
    
    %%%Calc Error Bar Lengths
    yneg = zp-mnz; ypos = mxz-zp;
    xneg = t-mnt(i)'; xpos = mxt(i)'-t; 
    xpos = xpos/2000*(max(rf_dist)); 
    xneg = xneg/2000*(max(rf_dist));
    
    %make matrix of results
    evt_list = [rf_dist(idx)' zp t]; 
  for ii = 1:size(DDA,1)
        if size(evt_list,1) >= 1;
           %errorbar(evt_list(:,1),evt_list(:,2),yneg,ypos,xneg,xpos,'k.','CapSize',10)
%            errorbar(evt_list(:,1),evt_list(:,2),yneg,ypos,'k.','CapSize',10)
           hold on
           scatter(evt_list(:,1),(evt_list(:,2)),60,evt_list(:,3),'d','filled','markeredgecolor','black');
           %errorbar(max(rf_dist)/10,(min(ylim)+ 20),50/2000*(max(rf_dist)),'horizontal','k')
           %text(max(rf_dist)/10,(min(ylim)+ 14),['50' char(176) ' C'], 'HorizontalAlignment','Center','FontSize',16);
        end
  end
%%%%%%%END Plot EVENTS%%%%%%%%%%%
