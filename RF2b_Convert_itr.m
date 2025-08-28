%Load radial receiver functions and convert to format CCP code can read this is
%only needed if you're importing SAC files.
clear all; close all
addpath ./Functions
sps = 10; dt = 1/sps; %new sampling rate
Stns = struct; Data = struct;% Create Structures to fill with station info and Data
Tcuts = [-5 85]; % time window to cut trace

files = dir('*.itr'); % find RF files
for i = 1:length(files)
    fname = [files(i).folder '/' files(i).name]; % find file name of individual RFs
    ddd = sac_read(fname); % read in sac file
    istn = cellstr(ddd.kstnm);	%get stn name
    intwrk = cellstr(ddd.knetwk); % get network name
    name = (char(istn)); % stn name as character array
    ntname = (char(intwrk));
    istn = [ntname '_' name];
    if isstrprop(name,'digit') % if name is a number add N to start of name
        istn = cellstr(['N' char(istn(i))]);
    end
    sr = ddd.delta; %get sample rate (Hz)
    rp = ddd.user(5); %get ray parameter
    tt = ddd.t; %get time vector of data
    rfd = ddd.d; %get RF values
    ntt = [Tcuts(1):1/sps:Tcuts(end)]; %create new time vector with new sample rate
    seis = interp1(tt,rfd,ntt); %interp to new sample rate and cut data
    if ~isfield(Stns,(char(istn))); % Assign Station info to structure if it doesn't exist already
        Stns.(char(istn)).Station_Data.Latitude = ddd.stla;
        Stns.(char(istn)).Station_Data.Longitude = ddd.stlo;  
        Stns.(char(istn)).Station_Data.Elevation = ddd.stel;
    end
    if rp ~= 0 %make sure ray parameter is assigned then put data into structure
        if ~isfield(Data,(char(istn)));
            j = 1
        else
            j = size(Data.(char(istn)).RFs_in_Time.RPs,2)+1;
        end
        Data.(char(istn)).RFs_in_Time.RPs(j) = rp.*(pi/180)./deg2km(1);
        Data.(char(istn)).RFs_in_Time.ID.RFs(j,:) = seis; 
        Data.(char(istn)).RFs_in_Time.ID.Time = ntt; 
        Data.(char(istn)).RFs_in_Time.BAZs(j) = ddd.baz;
    end
end

save('RF_converted.mat','Data','Stns')
  