function [eqd] = get_ph5(netwrk,stn,tstrt,tend,unamepwd)
base = './Functions/FetchData ';;
net = ['-N ' char(netwrk) ' '];
sta = ['-S ' char(stn) ' '];
tfrmt = 'yyyy-mm-ddTHH:MM:SS';
stime = ['-s ' datestr(tstrt,tfrmt) ' '];
etime = ['-e ' datestr(tend,tfrmt) ' '];
up = ['-a ' char(unamepwd(1)) ':' char(unamepwd(2)) ' '];
out = ['-o temp.mseed '];
sinfo = ['-m temp.meta '];
endinfo = '--timeseriesws http://service.iris.edu/ph5ws/dataselect/1/ --metadataws http://service.iris.edu/ph5ws/station/1/';
system([base net sta stime etime out sinfo up endinfo]); %download mseedfile

if exist('temp.meta') & exist('temp.mseed');
[X,I] = rdmseed('temp.mseed'); % read in mseed file
tend = dateshift(tend,'start','second');
tstrt = dateshift(tstrt,'start','second');

for i = 1:length(I);
    blks = I(i).XBlockIndex;
    XC = X(blks);
    
    %Get header data from data file
    snet = strtrim(XC(1).NetworkCode); % get network code
    snam = strtrim(X(1).StationIdentifierCode); % get stn code
    component = XC(1).ChannelIdentifier; % get component code
    sample_rate = unique([XC.SampleRate]); %get sample rate
    
    %Read from metadata file to get that info
    meta = readcell('temp.meta','Delimiter','|','DatetimeType','text','FileType','text');
    rw = find(strcmp(meta(:,4),component)); % fine component row
    meta = meta([1 rw],:); % remove other info
      col = find(strcmp(meta(1,:),'lon')); % fine lon column
    longitude = cell2mat(meta(2,col));
      col = find(strcmp(meta(1,:),'lat')); % fine lat column
    latitude = cell2mat(meta(2,col));   
      col = find(strcmp(meta(1,:),'depth')); % fine lat column
    depth = cell2mat(meta(2,col));
      col = find(strcmp(meta(1,:),'elev')); % fine lat column
    elevation = cell2mat(meta(2,col));
      col = find(strcmp(meta(1,:),'azimuth')); % fine lat column
    azimuth = cell2mat(meta(2,col));    
      col = find(strcmp(meta(1,:),'dip')); % fine lat column
    dip = cell2mat(meta(2,col)); 
      col = find(strcmp(meta(1,:),'instrument')); % fine lat column
    instrument = cell2mat(meta(2,col));
      col = find(strcmp(meta(1,:),'scale')); % fine lat column
    scale = cell2mat(meta(2,col));   
      col = find(strcmp(meta(1,:),'scalefreq')); % fine lat column
    scalefreq = cell2mat(meta(2,col));   
      col = find(strcmp(meta(1,:),'scaleunits')); % fine lat column
    scaleunits = cell2mat(meta(2,col));

    
    % Read Data from temp.mseed
    st = [cat(1,XC.t)]; % get time vector for cmpnt
    sd = [cat(1,XC.d)]; % get data vector for cmpnt
    [st,ind] = sort(st); % sort due to find extra data arriving
    sd = sd(ind);
    otime = st;
    st = datetime(st.*86400,'convertfrom','epochtime','Epoch','-0001-12-31 00:00:01','TicksPerSecond',1);% convert to time format data
    ind = find(st <= tend & st >= tstrt); % find extra data
    sd  = sd(ind); % cut out extra data
    st = st(ind);% cut out extra data
    otime = otime(ind); 
    stime = otime(1); 
    etime = otime(end);
    sample_count = length(st); %get sample count

   
    %Pole and zeros for fairfield 5hz
    if(contains(instrument,'Fairfield') & contains(instrument,'5Hz'))
        pz.k = [scale.*2.*pi.*scalefreq];
        pz.z = [0; 0];
        pz.p = [-2.19911E+01  -2.24355E+01i;
        -2.19911E+01  +2.24355E+01i];
    end
    trace.data = sd; trace.time = st; trace.sampleRate = sample_rate; % create trace file to remove instrument response
    sf = Remove_RESP(trace,.25,6.5,pz);  
    eqd(i).network = snet;
    eqd(i).station = stn;
    eqd(i).latitude = latitude;
    eqd(i).longitude = longitude;
    eqd(i).elevation = elevation;
    eqd(i).depth = depth;
    eqd(i).dip = dip;
    eqd(i).azimuth = azimuth;
    eqd(i).sampleCount = sample_count;
    eqd(i).sampleRate = sample_rate;
    eqd(i).channel = component;
    eqd(i).startTime = stime;
    eqd(i).endTime = etime;
    eqd(i).stime = st(1);
    eqd(i).etime = st(end);
    eqd(i).odata = -sd; %- accounts for reversed nodal polarity relative to BB
    eqd(i).data = -sf;% - accounts for reversed nodal polarity relative to BB
end
system(['rm temp.mseed temp.meta']);
else eqd = [];
end
end
    