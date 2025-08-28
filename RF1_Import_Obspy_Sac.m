close all; clear all;
addpath ./Functions

%Read off text file downloaded from https://ds.iris.edu/gmap/
    file = '6H.txt'
    data_directory = './RAW_Data'
    Squery = readtable(file,'CommentStyle','#');

    j = 1;
for i = 1:size(Squery,1);
    net = char(Squery{i,1});
    if isnumeric(Squery{i,2})
        stn = num2str(Squery{i,2});
    else
        stn = char(Squery{i,2});
    end
    if isfolder([data_directory '/' net '.' stn '/'])
        if ~isstrprop(net,'digit') %test if stn name is number, if so add N to structure name
            Stns.(char([net '_' stn])).Station_Data.Latitude = Squery{i,3};
            Stns.(char([net '_' stn])).Station_Data.Longitude = Squery{i,4};
            Stns.(char([net '_' stn])).Station_Data.Elevation = Squery{i,5};
        else 
            SC =['N' net '_' stn];
            Stns.(char(SC)).Station_Data.Latitude = Squery{i,3};
            Stns.(char(SC)).Station_Data.Longitude = Squery{i,4};
            Stns.(char(SC)).Station_Data.Elevation = Squery{i,5};
        end
    end
end
save('Stns.mat','Stns');


%%%Iterate through stations and find events in the time and distance window) %%%
slist = fields(Stns);
f1 = waitbar(0,'Importing Data');
for i = 1:length(slist)
    waitbar((i-1)/length(slist),f1,['Importing Station ' num2str(i) ' of ' num2str(length(slist))]);
        sinfo = strsplit(string(slist(i)),'_');
        net = char(sinfo(1));
        stn = char(sinfo(2));
        %flist = dir(([data_directory '/' net '.' stn '/' stn '*']));
        flist = dir(([data_directory '/' net(2:end) '.' stn '/' stn '*']));
        uflist = vertcat(flist.name);
        flist = unique(string(uflist(:,1:end-4)));
        for j = 1:length(flist)
            sfiles = dir([data_directory '/' net(2:end) '.' stn '/' char(flist(j)) '*']);
            if length(sfiles) ~= 3
                continue
            end
            for k = 1:3
                tr = sac_read([sfiles(k).folder '/' sfiles(k).name]);
                mday = datetime(2018,0,0)+days(tr.nzjday);
                eqd(k).network = strtrim(tr.kinst);
                eqd(k).station = stn;
                eqd(k).location = [];
                eqd(k).channel = strtrim(tr.kcmpnm);
                eqd(k).quality = [];
                eqd(k).latitude = tr.stla;
                eqd(k).longitude = tr.stlo;
                eqd(k).elevation = tr.evel;
                eqd(k).depth = [];
                eqd(k).azimuth = tr.cmpaz;
                eqd(k).dip = tr.cmpinc;
                eqd(k).sensitivity = [];
                eqd(k).sensitivityFrequency = [];
                eqd(k).instrument = [];
                eqd(k).sensitivityUnits = [];
                eqd(k).data = tr.d;
                eqd(k).sampleCount = length(tr.d);
                eqd(k).sampleRate = 1./(tr.delta);
                eqd(k).stime = datetime(tr.nzyear,month(mday),day(mday),tr.nzhour,tr.nzmin,tr.nzsec,tr.nzmsec); 
                eqd(k).etime = eqd(k).stime+seconds(tr.E);
                eqd(k).startTime = days(eqd(k).stime - datetime(0,1,1));
                eqd(k).endTime = days(eqd(k).etime - datetime(0,1,1));
                eqd(k).sacpz = [];
                eqd(k).eqinfo = tr.kevnm;
                eqd(k).PRP = tr.user(3);
                eqd(k).SRP = [];
                eqd(k).BAZ = tr.baz;
            end
            id = strtrim(tr.kevnm);
            EQ_Data.(id) = eqd;
            clear eqd id;
        end
        
    if exist('EQ_Data','var')
      save(['./Station_Data/' net stn '_RF_Data.mat'],'EQ_Data')
      clear EQ_Data
    end
end
close(f1);



                
  
            
        
        
