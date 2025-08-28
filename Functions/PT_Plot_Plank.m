function [] = PT_Plot(ctk,x,y,z,xd,yd,zd,rf_dist,file,binsz)
    pt = load(file);
   
    xpt = xd(1,:)'; ypt = yd(1,:)'; % define x and y pts
    
    
    evt_list = [];
    [idx d] = knnsearch([xpt ypt],[pt(:,2) pt(:,1)]); % find distances between samples and transect
    i = find(d < binsz/111);
    idx = idx(i); %choose samples within range
    d = d(i); % choose only distances within range
    pres = nan(size(z));
    pres(z >= -ctk) = (2700*-z(z >= -ctk)*1000*9.81);
    pres(z < -ctk) = (3300*((-z(z < -ctk))-ctk)*1000*9.81)+ctk*2700*1000*9.81;
    zp = zd(:,idx);
    pres = pres';
    zp = zd(:,idx);
    
    %calculate pressure in places where not previously constrained
    zp = (-zp)-min(-zp); 
    pres(isnan(pres)) =  zp(isnan(pres)).*1000*2700*9.81;% Pressure in crust
    
    %find closest pressure to PT estimates
    pi = find(abs(pres-pt(i,4)'*10^9) == min(abs(pres-pt(i,4)'*10^9)));
    zp = zd(pi); % depth
    t = pt(i,3); %temp
    
    
    %make matrix of results
    evt_list = [rf_dist(idx)' zp t]; 
        if size(evt_list,1) >= 1;
           hold on
           scatter(evt_list(:,1),(evt_list(:,2)),90,evt_list(:,3),'p','filled','markeredgecolor','black');
        end
%%%%%%%END Plot EVENTS%%%%%%%%%%%
