function [CCP HITS BSTD] = Grid_data(xpierce,ypierce,seis,x,y,z,bin_sz,bin_sm,gauss_width,bstrp);
 [X,Y] =  meshgrid(x,y); x = reshape(X,[numel(X),1]); y = reshape(Y,[numel(Y),1]); % create list of x y points
  CCP = NaN(length(z),length(x)); HITS = NaN(length(z),length(x)); BSTD = NaN(length(z),length(x)); 
  f = waitbar(0,'Building CCP Grid');
  for i=1:length(z); 
    waitbar(i/length(z),f,'Building CCP Grid') 
     [nrby,D] = (rangesearch([xpierce(i,:);ypierce(i,:)]',[x y],bin_sz*bin_sm));
      parfor j = 1:length(x)
      if (not(cellfun('isempty',nrby(j))))
        ind = cell2mat(nrby(j));
        if size(ind,2) > 0;
           DD = cell2mat(D(j));  % setup to weight RFs by distance from center of bin
           gaus = normpdf(DD,0,gauss_width);
           seisd = squeeze(seis(i,ind));
           gaus(isnan(seisd)) =[];
           seisd(isnan(seisd)) =[]; 
           if length(seisd ) < 5
              continue
           end
           [bootstat,~] = bootstrp(bstrp,@wmean,seisd,gaus);%,Options=statset(UseParallel=true));
           %pd = fitdist((squeeze(seis(i,ind))'),'normal');
           %pd = fitdist(seisd','Kernel','Kernel','normal');
           %pd_eval = linspace(-1,1,1000);
           %kd = kde(seisd,EvaluationPoints=pd,Kernel="box")
           %pind = find(pd.pdf(pd_eval) == max(pd.pdf(pd_eval)));
           %CCP(i,j) = sum(seisd.*gaus)./sum(gaus);
           % ci95 = paramci(pd);
           %BSTD(i,j) = std(seisd);
           %CCP(i,j) = pd.mean;
           %CCP(i,j) = trimmean(seisd,20);
           BSTD(i,j) = std(bootstat);
           %BSTD(i,j) = std(seisd);
           CCP(i,j) = mean(bootstat);%median(seis(i,ind));
           HITS(i,j) = length(seisd);
        else
          BSTD(i,j) = NaN;
          CCP(i,j) = NaN;
          HITS(i,j) = NaN;
        end
      else
        CCP(i,j) = NaN;
        HITS(i,j) = NaN;
        BSTD(i,j) = NaN;
      end;
    end;
 end;
CCP = reshape(CCP,[length(z),size(X)]);
HITS = reshape(HITS,[length(z),size(X)]);
BSTD = reshape(BSTD,[length(z),size(X)]);
close(f);
end

function tmean = wmean(sample,weight)
     oind = ~isoutlier(sample,"median");
     if isempty(sample(oind))
         return
     end
     tmean = sum(sample(oind).*weight(oind))./sum(weight(oind));
end