function [CCP HITS BSTD GSIZE] = Grid_data(xpierce,ypierce,seis,x,y,z,bin_min,bin_max,min_traces,gauss_width,bstrp,tnames);
 [X,Y] =  meshgrid(x,y); x = reshape(X,[numel(X),1]); y = reshape(Y,[numel(Y),1]); % create list of x y points
  CCP = NaN(length(z),length(x)); HITS = NaN(length(z),length(x)); BSTD = NaN(length(z),length(x)); GSIZE = NaN(length(z),length(x));
  f = waitbar(0,'Building CCP Grid');
  for i=1:length(z); 
    waitbar(i/length(z),f,'Building CCP Grid')
     [nrby,D] = (rangesearch([xpierce(i,:);ypierce(i,:)]',[x y],bin_max));
      parfor j = 1:length(x)
        if (not(cellfun('isempty',nrby(j))))
            ind = cell2mat(nrby(j)); 
            if size(ind,2) < min_traces;
                continue
            end
            DD = cell2mat(D(j));  % setup to weight RFs by distance from center of bin
            seisd = squeeze(seis(i,ind));
            anamesd = (tnames(ind));
            DD(isnan(seisd)) =[];
            seisd(isnan(seisd)) = [];
            if isempty(seisd)
                continue
            end
            [seis_s, gaus_s, bsearch] = stn_average(seisd,DD,anamesd,min_traces,3,bin_min,bin_max,gauss_width);
            % [seis_s, gaus_s, bsearch] = trace_average(seisd,DD,min_traces,bin_min,bin_max,gauss_width);
            if isnan(seis_s)
                continue
            end
           [bootstat,~] = bootstrp(bstrp,@wmean,seis_s,gaus_s);%;,Options=statset(UseParallel=true));
           %pd = fitdist((squeeze(seis(i,ind))'),'normal');
           %pd = fitdist(seisd','Kernel','Kernel','normal');
           %pd_eval = linspace(-1,1,1000);
           %kd = kde(seisd,EvaluationPoints=pd,Kernel="box")
           %pind = find(pd.pdf(pd_eval) == max(pd.pdf(pd_eval)));
           %CCP= pd_eval(pind);
           %CCP(i,j) = sum(seisd.*gaus)./sum(gaus);
           % ci95 = paramci(pd);
           %BSTD(i,j) = std(seisd);
           % CCP(i,j) = pd.mean;
           %CCP(i,j) = mean(seisd,'omitnan');
           BSTD(i,j) = std(bootstat);
           %BSTD(i,j) = std(seisd);
           CCP(i,j) = mean(bootstat);%median(seis(i,ind));
           HITS(i,j) = length(seisd);
           GSIZE(i,j) = bsearch;
        else
          BSTD(i,j) = NaN;
          CCP(i,j) = NaN;
          HITS(i,j) = NaN;
          GSIZE(i,j) = NaN;
        end
    end;
 end;
CCP = reshape(CCP,[length(z),size(X)]);
HITS = reshape(HITS,[length(z),size(X)]);
BSTD = reshape(BSTD,[length(z),size(X)]);
GSIZE = reshape(GSIZE,[length(z),size(X)]);
close(f);
end

function tmean = wmean(sample,weight)
     tmean = sum(sample.*weight)./sum(weight);
end

function [seis_s, gaus_s, bsearch] = trace_average(seisd,DD,min_traces,bin_min,bin_max,gauss_width)
    bsearch = bin_min;
    bind = find(DD < bsearch);
    while ~all([length(bind) > min_traces]) & bsearch < bin_max;
        bind = find(DD < bsearch);
        bsearch = bsearch + 1;
    end
    DD = DD(bind);
    seis_s = seisd(bind);
    gaus_s = normpdf(DD,0,gauss_width);
    if length(seis_s) < 2
      seis_s = nan;
      return
    end
end


function [seis_s, gaus_s, bsearch] = stn_average(seisd,DD,anamesd,min_traces,min_sta,bin_min,bin_max,gauss_width)
stblen = []; stbins = []; namesd = [];
    bsearch = bin_min;
    bind = find(DD < bsearch);
    while ~all([length(bind) > min_traces stblen > min_sta]) & bsearch < bin_max;
        bind = find(DD < bsearch);
        namesd = anamesd(bind);
        stbins = unique(namesd);
        stblen = length(stbins);
        bsearch = bsearch + 1;
    end
    namesd = anamesd(bind);
    stbins = unique(namesd);
    stblen = length(stbins);
    DD = DD(bind);
    seisd = seisd(bind);
    gaus = normpdf(DD,0,gauss_width);
    seis_s = [];
    gaus_s = [];
    nitr = 1;
    %tnum = ceil(length(bind)/stblen);
    tnum = [];
    for stinds = 1:stblen
        stind = find(strcmp(stbins(stinds),namesd));
        tnum = [tnum length(stind)];
    end
    tnum = ceil(median(tnum));
    for stinds = 1:stblen
        stind = find(strcmp(stbins(stinds),namesd)); 
        if (length(seisd(stind))) < tnum
            seis_s = [seis_s seisd(stind)];
            gaus_s = [gaus_s gaus(stind)];
        else
            svec = double(seisd(stind));
            sids = randperm(length(svec),tnum);
            gvec = gaus(stind);
            seis_s = [seis_s svec(sids)];
            gaus_s = [gaus_s gvec(sids)];
        end
         % if std(seisd(stind)) > .5 | length(seisd(stind)) < 2;
         %      continue
         %end 
         %seis_s(nitr) = median(seisd(stind));
         %gaus_s(nitr) = median(gaus(stind));
         %nitr = nitr + 1;
    end
      if length(seis_s) < 2
          seis_s = nan;
          return
     end
end