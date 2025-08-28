function [stn_list,sdist] = Stn_Line(x_val,y_val,rf_dist,Stns,binsz)
    slist = fields(Stns);
    for i = 1:length(slist)
        stn(i,:) = [Stns.(char(slist(i))).Station_Data.Longitude Stns.(char(slist(i))).Station_Data.Latitude  Stns.(char(slist(i))).Station_Data.Elevation]; 
    end
    evt_list = [];
    x = x_val'; y = y_val';
    [idx d] = knnsearch([x y],[stn(:,1) stn(:,2)]);
    i = find(d < binsz/111);
    idx = idx(i); d = d(i);
    stn_list = slist(i);
    sdist = [rf_dist(idx)'];
end