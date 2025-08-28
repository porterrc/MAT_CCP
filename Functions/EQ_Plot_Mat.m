function [] = EQ_Plot(x_val,y_val,rf_dist,binsz,file)
    EQ_Data = load(file);
    eq = [[EQ_Data.Event_Data.PreferredLongitude]'  [EQ_Data.Event_Data.PreferredLatitude]'  -[EQ_Data.Event_Data.PreferredDepth]']; 
    %x_err = mean(EQ_Data(:,5:6),2); z_err = (EQ_Data(:,7));
    evt_list = [];
    x = x_val'; y = y_val';
    [idx d] = knnsearch([x y],[eq(:,1) eq(:,2)]);
    i = find(d < binsz/111);
    idx = idx(i); d = d(i); z = eq(i,3);
    evt_list = [x(idx) y(idx) z];
    %err_lst = [x_err(i) z_err(i)];
    evt_list = [rf_dist(idx)' z];
    if size(evt_list,1) >= 1;
         plot(evt_list(:,1),(evt_list(:,2)),'o','MarkerFaceColor',[0.3010 0.7450 0.9330],'MarkerSize',20,'MarkerEdgeColor','k');
         %errorbar(evt_list(:,1),evt_list(:,2),err_lst(:,2),err_lst(:,2),err_lst(:,1),err_lst(:,1),'ko','MarkerFaceColor','b','MarkerSize',6)
         
    end
%%%%%%%END Plot EVENTS%%%%%%%%%%%
