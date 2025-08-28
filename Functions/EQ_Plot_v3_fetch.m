function [eq] = EQ_Plot(Start_Time,End_Time,Min_Mag,LL,cat,x_val,y_val,rf_dist,binsz)
    %Inputs are Start Time, End Time, Minimum Magitude, [Lat Limts Lon
    %Limits, catalog (ISC,NEIC%PDE, or GCMT).
    disp('Fetching dem EQs....Ohhh Yeeeahhh') 
    Event_Data = irisFetch.Events('StartTime' , Start_Time, 'EndTime', End_Time, 'minmag', Min_Mag,'boxcoordinates',LL,'catalog',cat);
    eq = [ cell2mat({Event_Data.PreferredLongitude})'  cell2mat({Event_Data.PreferredLatitude})'  -cell2mat({Event_Data.PreferredDepth})']; 
    evt_list = [];
    x = x_val; y = y_val;
    [idx d] = knnsearch([x y],[eq(:,1) eq(:,2)]);
    i = find(d < binsz/111);
    idx = idx(i); d = d(i); z = eq(i,3);
    evt_list = [x(idx) y(idx) z];
    evt_list = [rf_dist(idx)' z];
    if size(evt_list,1) >= 1;
         plot(evt_list(:,1),(evt_list(:,2)),'ko','MarkerFaceColor','w','MarkerSize',6);
         %errorbar(evt_list(:,1),evt_list(:,2),err_lst(:,2),err_lst(:,2),err_lst(:,1),err_lst(:,1),'ko','MarkerFaceColor','b','MarkerSize',6)
         
    end
%%%%%%%END Plot EVENTS%%%%%%%%%%%
