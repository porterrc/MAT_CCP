function [] = Ray_Path_Plot(x_val,y_val,rf_dist,binsz,xpierce,ypierce,z)
    rp_plot = nan(size(xpierce));   
    for i = 1:length(z)
        rp = [xpierce(i,:)'  ypierce(i,:)']; 
        evt_list = [];
        x = x_val'; y = y_val';
        [idx d] = knnsearch([x y],[rp(:,1) rp(:,2)]);
        ind = find(d >= binsz/111);
        ind0  = find(isnan(d));
        %idx(ind) = nan;
        rp_plot(i,:) = rf_dist(idx);
        rp_plot(i,ind) = nan;
        rp_plot(i,ind0) = nan;
    end
         plot(rp_plot,z'.*ones(size(rp_plot)),'k')
    end
%%%%%%%END Plot EVENTS%%%%%%%%%%%