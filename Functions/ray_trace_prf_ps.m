%%%need to work on veloictdy model import 3/17/23
function [z, xpierce,ypierce,baz, seisout, sd_name]  = ray_trace(Data,sname,stnx,stny,stel,dzi,z_max,sx,sy,sz,svmodel,pvmodel,vmodel,minv,mxamp,mnamp,blim)
 decon = 'ID'
 ref_angle = 15; %reference incident angle for amplitude scaling
 vel = dlmread(vmodel,'', 2, 0); % read in velocity mode for areas outside shear model
 z = (-round(max(stel)*1/dzi)*dzi):dzi:z_max; % creatw depth vector
 tmp = vel(:,1);
 [~,ind,~] = unique(tmp); tmp(ind) = -1; ind = find(tmp > 0); % account for duplicate values
 vel(ind,1) = vel(ind,1)+0.0001; % account for duplicate values
 vpi = interp1(vel(:,1),vel(:,2),z(1:end)+0.5*dzi); % interpolate vp and vs to depth
 vsi = interp1(vel(:,1),vel(:,3),z(1:end)+0.5*dzi); % interpolate vp and vs to depth 
 [SX SY] = meshgrid(sx,sy); %convert shear velocity locations to vector
 SX = reshape(SX,[numel(SX) 1]); SY = reshape(SY,[numel(SX) 1]);
 svect = reshape(svmodel,[length(SX),size(svmodel,3)]);
 pvect = reshape(pvmodel,[length(SX),size(svmodel,3)]);
 k = 1;
for i = 1:length(sname);
    fd_name = char(sname(i));
    % if contains(fd_name,'_')
    %     finfo = strsplit(fd_name,'_');
    %     fd_name = ['N' char(finfo(2))];
    % end
    if isfield(Data,fd_name);      
        p = Data.(fd_name).RFs_in_Time.RPs;
        backaz = Data.(fd_name).RFs_in_Time.BAZs;
        if ~isfield(Data.(fd_name).RFs_in_Time,decon)
            continue
        end
        seis = Data.(fd_name).RFs_in_Time.(decon).RFs;
        time = Data.(fd_name).RFs_in_Time.(decon).Time;
        stx = stnx(i); sty = stny(i); % station locations on grid
        se = round(stel(i)./dzi)*dzi; % rounded station location to dzi increment
        [s_ind d]= knnsearch([SX,SY],[stx,sty]);
        if d <= 50 % only use if svelocity within 50 km of sample point
            vsu = svect(s_ind,:)'; %
            vsu(vsu <= minv) = min(vsu(vsu>minv));
            vs = interp1(sz,vsu,z(1:end)+0.5*dzi,'linear','extrap');
            vs(vs < min(vsu)) = min(vsu);
            vs = fillmissing(vs,'nearest');
            vpu = pvect(s_ind,:)'; %
            vpu(vpu <= minv) = min(vpu(vpu>minv));
            vp = interp1(sz,vpu,z(1:end)+0.5*dzi,'linear','extrap');
            vp(vp < min(vpu)) = min(vpu);
            vp = fillmissing(vp,'nearest');
            %vp = vs*1.75; % set vp/vs to 1.75. this could also be queried
        else
            vs = vsi;
            vs= fillmissing(vs,'nearest');
            vp = vpi;
            vp= fillmissing(vp,'nearest');
        end
        if min(isnan(vs)) == 1;
            vs = vsi;
            vs= fillmissing(vs,'nearest');
            vp = vpi;
            vp= fillmissing(vp,'nearest');
        end
            
    %    for j = find(Data.(char(sname(i))).RFs_in_Time.(decon).Grade <= 1);
        for j = 1:size(Data.(fd_name).RFs_in_Time.(decon).RFs,1);
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
              if max(seis(j,:)) > mxamp || min(seis(j,:)) < mnamp;
                continue
              end
                if backaz(j) > blim(1) &&  backaz(j) < blim(2);
                  eci = (((-min(z)-se)/dzi+1)); % elevation correction indices            
                  seist = interp1(time,seis(j,:),tmid,'linear'); % calulate the amplitude at the middle of layers
                  seist = [NaN(1,eci) seist]; seist = seist(1:end-eci); % add and remove values to acct for elevation
                  seis_amp = seist.*a_scale;
                  seisout(:,k) = seist.*a_scale; % mulitply by incident angle scalar to correct to ref angle
                  qb = qb(ind); qa = qa(ind);
                  dh = NaN(size((qb))); hpos(:,k) = NaN(size((qb))); % Define Variables
                  dh(eci:end) = dzi*p(j)./qb(eci:end); % horizontal offset in each depth interval
                  hpos((eci:end),k) = [0, cumsum(dh(eci:end-1))]; % total horizontal offset at depth
                  ew(:,k) = hpos(:,k).*sind(backaz(j)); % ew offset
                  ns(:,k) = hpos(:,k).*cosd(backaz(j)); % ns offsest
                  baz(k) = backaz(j);
                  xpierce(:,k) = ew(:,k) + stx; 
                  ypierce(:,k) = ns(:,k) + sty;
                  sd_name(k) = {fd_name};
                  k = k + 1;
                end
        end
    end
end
xpierce(isinf(xpierce)) = NaN;
ypierce(isinf(ypierce)) = NaN;
