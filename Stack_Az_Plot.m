%stack_azz_plot.m
%Plots back-azimuths for radial and tangential receiver functions.
clear all; close all


addpath('./Functions')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%USER CONTROLLED INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ifile = 'RF_5_QC_Merge.mat';      %filename start
idir  = 'azz_plots/';   %output directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mkdir(idir);            %make directory
dt = .1
fprintf('==========Load Dataset========== \n')
load(ifile)
%window seismograms
%binning parameters
bins = 0:5:360; nbins = size(bins,2); brng = 5;

%%%Get Station Names%%%
snames = fields(Data);
ns = length(snames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%iterate through stations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for s = 1:ns
   figure(1); clf
   for i = [1 2]
      subplot(1,2,i)
      tt = Data.(char(snames(s))).RFs_in_Time.ID.Time';
      npts = length(tt);
      ntra = length(Data.(char(snames(s))).RFs_in_Time.RPs);
      binstk = NaN(npts,nbins);
      if i == 1; 
          sss = -Data.(char(snames(s))).RFs_in_Time.ID.RFs';
            mind = find(max(abs(sss),[],1) > 1.5);
            if size(mind,2) ~=0;
                sss(:,mind) = [];
                Data.(char(snames(s))).RFs_in_Time.BAZs(mind) = [];
                Data.(char(snames(s))).RFs_in_Time.ID.RFts(mind,:) = [];
            end
      end
      if i == 2; 
          sss = -Data.(char(snames(s))).RFs_in_Time.ID.RFts'; 
        mind = find(max(abs(sss),[],1) > 2.5);
        if size(mind,2) ~=0
            sss(:,mind) = [];
            Data.(char(snames(s))).RFs_in_Time.BAZs(mind) = [];
        end
      end
      for j = 1:nbins
         wnt = find(Data.(char(snames(s))).RFs_in_Time.BAZs > bins(j)-brng & Data.(char(snames(s))).RFs_in_Time.BAZs < bins(j)+brng);
         if size(wnt,2) > 1
            binstk(:,j) = median(sss(:,wnt)')'; binstk(npts,j) = 0;
         elseif size(wnt,2) == 1
            binstk(:,j) = sss(:,wnt); binstk(npts,j) = 0;
         end
         bb = (binstk(:,j)*20)+bins(j);
         plot(tt,bb,'k'); hold on
         ff = find(bb <= bins(j)); tbb = bb; tbb(ff) = bins(j);
         fill([tt' min(tt) min(tt)],[tbb' tbb(1) tbb(1)],'r')
         ff = find(bb >= bins(j)); bbb = bb; bbb(ff) = bins(j); 
         fill([tt' min(tt) min(tt)],[bbb' bbb(1) bbb(1)],'b')
         clear wnt
      end
      xlabel('time after P-wave (secs)'); ylabel('back-azimuth (degrees)')
      axis([-1 10 -15 375])
      set(gca,'XTick',0:2:36); set(gca,'YTick',0:40:360)
      if i == 1; title(['P and PP data ' snames{s} ' ' num2str(ntra) ' radial traces'], 'Interpreter', 'none'); end
      if i == 2; title(['P and PP data ' snames{s} ' ' num2str(ntra) ' tangential traces'], 'Interpreter', 'none'); end
      grid on
   end
   set(gcf,'paperorientation','portrait','paperposition',[0 0 8.5 11])
   eval(['print -dpng ' idir 'stack_azz_a' snames{s} '.png']);
%  pause
%   eval(['print -djpeg100 ' idir 'stack_azz_a' num2str(nbp(1)) '_' char(stns(s))]);
end
%quit
