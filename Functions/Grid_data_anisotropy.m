function [CCP CCP_RAW BAZ_RAW HITS BSTD] = Grid_data_v4(xpierce,ypierce,seis,baz,x,y,z,bin_sz,bin_sm,gauss_width,bstrp);
 [X,Y] =  meshgrid(x,y); x = reshape(X,[numel(X),1]); y = reshape(Y,[numel(Y),1]); % create list of x y points
  CCP = NaN(length(z),length(x)); HITS = NaN(length(z),length(x)); BSTD = NaN(length(z),length(x)); 
  CCP_RAW = NaN(length(z),length(x),size(seis,2)); BAZ_RAW = NaN(length(z),length(x),size(seis,2));
  f = waitbar(0,'Building CCP Grid');
  for i=1:length(z); 
    waitbar(i/length(z),f,'Building CCP Grid') 
     [nrby,D] = (rangesearch([xpierce(i,:);ypierce(i,:)]',[x y],bin_sz*bin_sm));
      for j = 1:length(x)
      if (not(cellfun('isempty',nrby(j))))
        ind = cell2mat(nrby(j));
        if size(ind,2) > 20;
          DD = cell2mat(D(j));  % setup to weight RFs by distance from center of bin
          gaus = normpdf(DD,0,gauss_width);
          brng = 0:10:360;
          sbaz = baz(ind);
          bseis = (seis(i,ind));
          % bwseis = nan;
          % bb = 1;
          % for b = 1:length(brng)-1
          %       bind = find(sbaz >= brng(b) & sbaz <= brng(b+1));
          %       if length(bind) == 0;
          %           continue
          %       end
          %       bnseis= bseis(bind);
          %       bngaus = gaus(bind);
          %       nseis = ~isnan(bnseis);
          %           if nseis == 0
          %       end
          %       bwseis(bb) = mean(bnseis(nseis).*bngaus(nseis))./sum(bngaus(nseis));
          %       bb = bb + 1;
          %   end
          %   if length(bwseis) == 1
          %       continue
          %   end
            [bootstat,~] = bootstrp(bstrp,@mean,bseis);
            CCP_RAW(i,j,1:length(bseis)) = bseis;
            BAZ_RAW(i,j,1:length(sbaz)) = sbaz;
            BSTD(i,j) = std(bootstat);
            CCP(i,j) = mean(bootstat);%median(seis(i,ind));
            HITS(i,j) = length(seis(i,ind));
            clear bwseis
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
cind = squeeze(all(isnan(CCP_RAW),[1 2]));
CCP_RAW(:,:,cind) = [];
BAZ_RAW(:,:,cind) = [];
CCP_RAW = reshape(CCP_RAW,[length(z),size(X),size(BAZ_RAW,3)]);
BAZ_RAW = reshape(BAZ_RAW,[length(z),size(X),size(BAZ_RAW,3)]);
close(f);
