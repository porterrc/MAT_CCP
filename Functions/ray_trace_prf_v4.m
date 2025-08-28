function [z, xpierce,ypierce, baz, seisout]  = ray_trace(Data,sname,stnx,stny,stel,dzi,z_max,vmodel,mxamp,mnamp,blim)
 ref_angle = 20; %reference incident angle for amplitude scaling
 vel = dlmread(vmodel,'', 2, 0); % read in velocity mode
 z = (-round(max(stel)*1/dzi)*dzi):dzi:z_max; % creat depth vector
 tmp = vel(:,1);
 [~,ind,~] = unique(tmp); tmp(ind) = -1; ind = find(tmp > 0); % account for duplicate values
 vel(ind,1) = vel(ind,1)+0.0001; % account for duplicate values
 vp = interp1(vel(:,1),vel(:,2),z(1:end)+0.5*dzi); % interpolate vp and vs to depth
 vs = interp1(vel(:,1),vel(:,3),z(1:end)+0.5*dzi); % interpolate vp and vs to depth 
 vp(isnan(vp)) = vel(1,2);
 vs(isnan(vs)) = vel(1,3); 
 k = 1;
for i = 1:length(sname);
    if isfield(Data,char(sname(i)));
        p = Data.(char(sname(i))).RFs_in_Time.RPs;
        backaz = Data.(char(sname(i))).RFs_in_Time.BAZs;
        seis = Data.(char(sname(i))).RFs_in_Time.ID.RFs;
        time = Data.(char(sname(i))).RFs_in_Time.ID.Time;
        sx = stnx(i); sy = stny(i); % station locations on grid
        se = round(stel(i)./dzi)*dzi; % rounded station location to dzi increment
    %    for j = find(Data.(char(sname(i))).RFs_in_Time.ID.Grade <= 1);
        for j = 1:size(Data.(char(sname(i))).RFs_in_Time.ID.RFs,1);
          qa = real(sqrt(1./vp.^2 - p(j).^2));
          qb = sqrt(1./vs.^2 - p(j).^2);
          i_ang = asind(p(j).*vs);
          a_scale = (ref_angle./i_ang); % scale to account for incident angle
          dt = dzi*(qb - qa); % time difference at each depth interval
          tbot = cumsum(dt); % get the time at bottom of each layer
          ttop = tbot - dt;  %get time at top of each layer
          tmid = mean([ttop; tbot],1);  % get time in middle of layer
          ind = find(tbot  <= max(time)); %get rid of layers deeper than trace "sees"
          tbot = tbot(ind); ttop= ttop(ind); tmid = tmid(ind); z = z(ind);
          if max(abs(seis(j,:))) < mxamp && min(seis(j,:)) > mnamp;
            if backaz(j) > blim(1) ||  backaz(j) > blim(2);
              eci = (((-min(z)-se)/dzi+1)); % elevation correction indices            
              seist = interp1(time,seis(j,:),tmid,'linear'); % calulate the amplitude at the middle of layers
              seist = [NaN(1,eci) seist]; seist = seist(1:end-eci); % add and remove values to acct for elevation
              seisout(:,k) = seist.*a_scale; % mulitply by incident angle scalar to correct to ref angle
              qb = qb(ind); qa = qa(ind);
              dh = NaN(size((qb))); hpos(:,k) = NaN(size((qb))); % Define Variables
              dh(eci:end) = dzi*p(j)./qb(eci:end); % horizontal offset in each depth interval
              hpos((eci:end),k) = [0, cumsum(dh(eci:end-1))]; % total horizontal offset at depth
              ew(:,k) = hpos(:,k).*sind(backaz(j)); % ew offset
              ns(:,k) = hpos(:,k).*cosd(backaz(j)); % ns offsest
              baz(k) = backaz(j);
              xpierce(:,k) = ew(:,k) + sx; 
              ypierce(:,k) = ns(:,k) + sy;
              k = k + 1;
            end
          end
        end
    end
end
xpierce(isinf(xpierce)) = NaN;
ypierce(isinf(ypierce)) = NaN;
