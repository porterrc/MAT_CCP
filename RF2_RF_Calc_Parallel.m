clear all; close all;
addpath ./Functions/

%%%%%Prep Data and Calculate RF%%%%%
files = ls('./Station_Data/*RF_Data.mat');
files = sort(strsplit(files,{'.mat','\t','\n'}));
files(1) = [];
rsr = 25; %resample rate
lp = 5; % lowpass
hp = .25; %highpass
Data = [];
for i = 1:length(files)
    baz = []; RPs = []; RFr = []; Err = []; RFt = []; Ert = []; Dsnr = []; data = [];
    %waitbar((i-2)/(length(files)-1),f1,['Processing Station ' num2str(i-1) ' of ' num2str(length(files)-1)]);
    SData = load([char(files(i)) '.mat']);
    EQ_Data = SData.EQ_Data;
    evs = fields(EQ_Data);
    if ~(isstrprop(EQ_Data.(char(evs(1)))(1).network(1),'digit'))
        s_name = [EQ_Data.(char(evs(1)))(1).network '_' EQ_Data.(char(evs(1)))(1).station];
    else
        s_name = ['N' EQ_Data.(char(evs(1)))(1).network '_' EQ_Data.(char(evs(1)))(1).station];
    end
        ev = 1;
        while isempty(s_name)
            ev = ev + 1;
            s_name = EQ_Data.(char(evs(ev)))(1).station;
        end
    disp([s_name ' ' char(datetime('now'))]);
    events = fields(EQ_Data);
    l = 1;
    for j= 1:length(events);
        %Sync, cut, and resample data
        stimes = cell2mat({EQ_Data.(char(events(j))).startTime})*86400;
        if(cellfun('isempty',{EQ_Data.(char(events(j))).station}) == 0)
        if (length({EQ_Data.(char(events(j))).station}) == 3)    
            etimes = cell2mat({EQ_Data.(char(events(j))).endTime})*86400;
            baz(l) = unique(cell2mat({EQ_Data.(char(events(j))).BAZ}));
            RPs(l) = unique(cell2mat({EQ_Data.(char(events(j))).PRP}));
            sr = cell2mat({EQ_Data.(char(events(j))).sampleRate});
            channels = string({EQ_Data.(char(events(j))).channel});
            new_st = max(stimes);
            new_et = min(etimes);
            ind = endsWith(channels,'Z');
            n=(numel([EQ_Data.(char(events(j)))(ind).data]));
            tsz = timeseries([EQ_Data.(char(events(j)))(ind).data(1:n)],[stimes(ind):1./sr(ind):stimes(ind)+(n-1)./sr(ind)]);
            ind = endsWith(channels,["E","2"]);
            n=(numel([EQ_Data.(char(events(j)))(ind).data]));
            tse = timeseries([EQ_Data.(char(events(j)))(ind).data(1:n)],[stimes(ind):1./sr(ind):stimes(ind)+(n-1)./sr(ind)]);
            tseaz = (([EQ_Data.(char(events(j)))(ind).azimuth]));
            ind = endsWith(channels,["N","1"]);
            n=(numel([EQ_Data.(char(events(j)))(ind).data]));
            tsn = timeseries([EQ_Data.(char(events(j)))(ind).data(1:n)],[stimes(ind):1./sr(ind):stimes(ind)+(n-1)./sr(ind)]);
            tsnaz = (([EQ_Data.(char(events(j)))(ind).azimuth]));
            %%cut data
            tsz = getsampleusingtime(tsz,stimes(1)+(50-10),stimes(1)+(50+100));
            tse = getsampleusingtime(tse,stimes(1)+(50-10),stimes(1)+(50+100));
            tsn = getsampleusingtime(tsn,stimes(1)+(50-10),stimes(1)+(50+100));
            if (any(tsz.data) & any(tsn.data) & any(tse.data));
                %%%sync data
                [tsz,tse] = synchronize(tsz,tse,'Uniform','Interval',1./rsr);
                [tse,tsn] = synchronize(tse,tsn,'Uniform','Interval',1./rsr);
                [tsn,tsz] = synchronize(tsn,tsz,'Uniform','Interval',1./rsr);
                chk = [length(tsz.data) length(tse.data) length(tsn.data)];
                while (length(unique(chk)) ~= 1)
                    [tsz,tse] = synchronize(tsz,tse,'Uniform','Interval',1./rsr);
                    [tse,tsn] = synchronize(tse,tsn,'Uniform','Interval',1./rsr);
                    [tsn,tsz] = synchronize(tsn,tsz,'Uniform','Interval',1./rsr);
                    chk = [length(tsz.data) length(tse.data) length(tsn.data)];
                end

                data(:,1) = tsz.Data; data(:,2) = tse.Data; data(:,3) = tsn.Data;
                L = size(data,1); % Build Taper
                Taper= tukeywin(L,0.05); % Build Taper
                data = (data.*Taper); % Taper trace at edges
                %%%rmean, rtrend, filter%%%
                data = data - mean(data,1);
                data = detrend(data);
                bpFilt = designfilt('bandpassiir','FilterOrder',4,'PassbandFrequency1',hp,'PassbandFrequency2',lp, 'SampleRate',rsr); % build  filter
                data = filter(bpFilt,data);
                data = (data.*Taper); % Taper trace at edges
                %%%%Rotate%%%%%
                [r t] = rotate_rt(data(:,3),data(:,2),baz(l),tsnaz) ; % n,e,baz
                %%%%Calculate RF%%%%%
                data(:,4) = r;
                data(:,5) = t;
                npow = mean(data(20.*rsr:30.*rsr,[1 4 5]).^2);
                spow = mean(data(40.*rsr:50.*rsr,[1 4 5]).^2);
                snr = mean(spow./npow);
                data(1:30*rsr,:) = [];
                L = size(data,1); % Build Taper
                Taper= tukeywin(L,0.05); % Build Taper
                data = (data.*Taper); % Taper trace at edges
                if size(data,1) > 1990;
                    %%%%%%Write out RF to File%%%%%
                    [itr rvar] = makeRFitdecon(data(:,4),data(:,1),1/rsr,0,30,10,2.5,1200,0.00001,0);
                    [itt tvar] = makeRFitdecon(data(:,5),data(:,1),1/rsr,0,30,10,2.5,1200,0.00001,0);
                    RFr(l,:) = itr; 
                    Err(:,l) = rvar(end);
                    RFt(l,:) = itt; 
                    Ert(:,l) = tvar(end);
                    Dsnr(l) = snr;
                    l = l + 1;
                end
            end
        end
        end
    data = [];
    end

    if ~isempty(RPs) 
    if isstrprop(s_name,'digit')
        s_name = ['N' s_name];
    end
    Data.(s_name).RFs_in_Time.RPs = RPs*pi/180./deg2km(1);
    Data.(s_name).RFs_in_Time.BAZs = baz;
    Data.(s_name).RFs_in_Time.ID.RFs = RFr;
    Data.(s_name).RFs_in_Time.ID.Error = Err;
    Data.(s_name).RFs_in_Time.ID.RFts = RFt;
    Data.(s_name).RFs_in_Time.ID.Errort = Ert;
    Data.(s_name).RFs_in_Time.ID.Time = [-10:1./rsr:30];
    Data.(s_name).RFs_in_Time.ID.SNR = Dsnr;
    baz = []; RPs = []; RFr = []; Err = []; RFt = []; Ert = []; Dsnr = []; 
    end
end
load Stns.mat
save('RF_mat.mat','Data','Stns')