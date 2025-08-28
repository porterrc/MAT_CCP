close all; clear all;
ifile = 'RF_5_Merge.mat'; % file name with RFs in it
ofile = 'RF_5_QC_Merge.mat'; % file name for output
mfit = 0.5 % maximum RF variance to accept
msnr = [1 1 0.5] % min snr to accept for v, r, and t components
amps = [-3 3] % min and max amplitudes to accept within a trace
iamp = 0 % minimum amplitude of first arrival
max_width = 4 % max size of a pulse in seconds
famax = 2; % max amplitude allowed greater than initial peak
%%%the next two lines use a peak finder to look for ringy traces
pkamp = 0.15; %min amp to use to look for pksand troughs
pkdist = .01; %min mean dist between pks and troughs
manual_qc = 1% set to one if you want to look at traces manually
load(ifile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Should write function to do this
fprintf('Pulling out Traces with High Variance \n')
    stns = fields(Data);
    for i = 1:length(stns)
        ind = find(Data.(char(stns(i))).RFs_in_Time.ID.Error > mfit);
        Data.(char(stns(i))).RFs_in_Time.RPs(ind) = [];
        Data.(char(stns(i))).RFs_in_Time.BAZs(ind) = [];
        %Data.(char(stns(i))).RFs_in_Time.SNR(ind,:) = [];
        Data.(char(stns(i))).RFs_in_Time.ID.RFs(ind,:) = [];
        Data.(char(stns(i))).RFs_in_Time.ID.Error(ind) = [];
        Data.(char(stns(i))).RFs_in_Time.ID.RFts(ind,:) = [];
        Data.(char(stns(i))).RFs_in_Time.ID.Errort(ind) = [];   
    end

% fprintf('Pulling out Traces with Low SNR \n')
%     stns = fields(Data);
%     for i = 1:length(stns)
%         if isempty(Data.(char(stns(i))).RFs_in_Time.SNR)
%             continue
%         end
%         ind = any((Data.(char(stns(i))).RFs_in_Time.SNR < msnr),2);
%         Data.(char(stns(i))).RFs_in_Time.RPs(ind) = [];
%         Data.(char(stns(i))).RFs_in_Time.BAZs(ind) = [];
%         Data.(char(stns(i))).RFs_in_Time.SNR(ind,:) = [];
%         Data.(char(stns(i))).RFs_in_Time.ID.RFs(ind,:) = [];
%         Data.(char(stns(i))).RFs_in_Time.ID.Error(ind) = [];
%         Data.(char(stns(i))).RFs_in_Time.ID.RFts(ind,:) = [];
%         Data.(char(stns(i))).RFs_in_Time.ID.Errort(ind) = [];   
%     end    
    
    
fprintf('Pulling out Traces with Amplitudes Outside Range \n')
    stns = fields(Data);
    for i = 1:length(stns)
        if isempty(Data.(char(stns(i))).RFs_in_Time.ID.RFs)
            continue
        end
        ind = find(min(Data.(char(stns(i))).RFs_in_Time.ID.RFs,[],2) < amps(1) | max(Data.(char(stns(i))).RFs_in_Time.ID.RFs,[],2) > amps(2));
        Data.(char(stns(i))).RFs_in_Time.RPs(ind) = [];
        Data.(char(stns(i))).RFs_in_Time.BAZs(ind) = [];
        %Data.(char(stns(i))).RFs_in_Time.SNR(ind,:) = [];
        Data.(char(stns(i))).RFs_in_Time.ID.RFs(ind,:) = [];
        Data.(char(stns(i))).RFs_in_Time.ID.Error(ind) = [];
        Data.(char(stns(i))).RFs_in_Time.ID.RFts(ind,:) = [];
        Data.(char(stns(i))).RFs_in_Time.ID.Errort(ind) = [];   
    end    
    
 fprintf('Pulling out Traces with First Arrivals Less than iamp\n')
    % stns = fields(Data);
    % for i = 1:length(stns)
    %     RF  = Data.(char(stns(i))).RFs_in_Time.ID.RFs;
    %     TV =  Data.(char(stns(i))).RFs_in_Time.ID.Time;
    %     if isempty(RF)
    %         continue
    %     end
    %     twin = find(TV > -.5 & TV < .5); % timewindow when 1st arrival expected 
    %     [~,ind] = max(abs(RF(:,twin)),[],2);
    %     ind = find(diag(RF(:,twin(ind))) < iamp);
    %     Data.(char(stns(i))).RFs_in_Time.RPs(ind) = [];
    %     Data.(char(stns(i))).RFs_in_Time.BAZs(ind) = [];
    %     %Data.(char(stns(i))).RFs_in_Time.SNR(ind,:) = [];
    %     Data.(char(stns(i))).RFs_in_Time.ID.RFs(ind,:) = [];
    %     Data.(char(stns(i))).RFs_in_Time.ID.Error(ind) = [];
    %     Data.(char(stns(i))).RFs_in_Time.ID.RFts(ind,:) = [];
    %     Data.(char(stns(i))).RFs_in_Time.ID.Errort(ind) = [];   
    % end
   

   fprintf('Pulling out Traces where First Arrival is not the Highest Amp \n')
    % stns = fields(Data);
    % for i = 1:length(stns)
    %     RF  = Data.(char(stns(i))).RFs_in_Time.ID.RFs;
    %     TV =  Data.(char(stns(i))).RFs_in_Time.ID.Time;
    %     if isempty(RF)
    %         continue
    %     end
    %     twin = find(TV > -.75 & TV < .75); % timewindow when 1st arrival expect
    %     [~,ind] = max((RF(:,twin)),[],2);
    %     mval = diag(RF(:,twin(ind)));
    %     tmax = max(abs(RF),[],2);
    %     ind = find((abs(mval) - abs(tmax)) < -famax);
    %     Data.(char(stns(i))).RFs_in_Time.RPs(ind) = [];
    %     Data.(char(stns(i))).RFs_in_Time.BAZs(ind) = [];
    %     %Data.(char(stns(i))).RFs_in_Time.SNR(ind,:) = [];
    %     Data.(char(stns(i))).RFs_in_Time.ID.RFs(ind,:) = [];
    %     Data.(char(stns(i))).RFs_in_Time.ID.Error(ind) = [];
    %     Data.(char(stns(i))).RFs_in_Time.ID.RFts(ind,:) = [];
    %     Data.(char(stns(i))).RFs_in_Time.ID.Errort(ind) = [];    
    % end  
    
fprintf('Pulling out Traces with Too Wide of Pulses \n')
    % stns = fields(Data);
    % for i = 1:length(stns)
    %     RF  = Data.(char(stns(i))).RFs_in_Time.ID.RFs;
    %     TV =  Data.(char(stns(i))).RFs_in_Time.ID.Time;
    %     if isempty(RF)
    %         continue
    %     end
    %     sr = TV(2)-TV(1);
    %     [b a] = find(diff(sign(RF')));
    %     ind = [];
    %     for av = 1:max(unique(a));
    %         aind = find(a == av);
    %         if(max(sr.*diff(b(aind)) > max_width));
    %             ind = [ind av];
    %         end
    %     end
    %     Data.(char(stns(i))).RFs_in_Time.RPs(ind) = [];
    %     Data.(char(stns(i))).RFs_in_Time.BAZs(ind) = [];
    %     %Data.(char(stns(i))).RFs_in_Time.SNR(ind,:) = [];
    %     Data.(char(stns(i))).RFs_in_Time.ID.RFs(ind,:) = [];
    %     Data.(char(stns(i))).RFs_in_Time.ID.Error(ind) = [];
    %     Data.(char(stns(i))).RFs_in_Time.ID.RFts(ind,:) = [];
    %     Data.(char(stns(i))).RFs_in_Time.ID.Errort(ind) = [];   
    %end    

    %uses peak finder to find mean distance between peaks and troughs
fprintf('Pulling Waves with lots of oscillations (ie ringy)\n')
    % stns = fields(Data);
    % for i = 1:length(stns)
    %     RF  = Data.(char(stns(i))).RFs_in_Time.ID.RFs;
    %     TV =  Data.(char(stns(i))).RFs_in_Time.ID.Time;
    %     if isempty(RF)
    %         continue
    %     end
    %     twin = find(TV > -.75 & TV < 30);
    %     if size(RF,1) > 0;
    %         ind = [];
    %         for j = 1:size(RF,1)
    %             [ppks,plocs] = findpeaks(RF(j,:),TV,'MinPeakHeight',.05);
    %             [npks,nlocs] = findpeaks(-RF(j,:),TV,'MinPeakHeight',.05);
    %             locs = mean(diff(sort([ plocs nlocs])));
    %             if locs < pkdist
    %                 ind = [ind j];
    %             end
    %         end
    %             Data.(char(stns(i))).RFs_in_Time.RPs(ind) = [];
    %             Data.(char(stns(i))).RFs_in_Time.BAZs(ind) = [];
    %            % Data.(char(stns(i))).RFs_in_Time.SNR(ind,:) = [];
    %             Data.(char(stns(i))).RFs_in_Time.ID.RFs(ind,:) = [];
    %             Data.(char(stns(i))).RFs_in_Time.ID.Error(ind) = [];
    %             Data.(char(stns(i))).RFs_in_Time.ID.RFts(ind,:) = [];
    %             Data.(char(stns(i))).RFs_in_Time.ID.Errort(ind) = [];  
    %         end
    % end
        
    
if manual_qc == 1;    
    fprintf('Manual QC \n')
    stns = fields(Data);
    for i = 1:length(stns)
        RF  = -Data.(char(stns(i))).RFs_in_Time.ID.RFs;
        TV =  Data.(char(stns(i))).RFs_in_Time.ID.Time;
      if size(RF,1) > 0;
        fig = figure
        plot(TV,RF)
        xlim([-5 20])
        
        title(char(stns(i)));
        datacursormode on  % Enable data cursor mode
        dcm_obj = datacursormode(fig);
        set(dcm_obj,'UpdateFcn',@myupdatefcn) % Set update function
        % Wait while the user to click
        disp('Click line to display a data tip, then press "Return"')
        pause 
        info_struct = getCursorInfo(dcm_obj); % Export cursor to workspace
        close(fig)
        while isfield(info_struct, 'Position')
           ind = find(RF(:,info_struct.DataIndex) == info_struct.Position(2));
            Data.(char(stns(i))).RFs_in_Time.RPs(ind) = [];
            Data.(char(stns(i))).RFs_in_Time.BAZs(ind) = [];
            %Data.(char(stns(i))).RFs_in_Time.SNR(ind,:) = [];
            Data.(char(stns(i))).RFs_in_Time.ID.RFs(ind,:) = [];
            Data.(char(stns(i))).RFs_in_Time.ID.Error(ind) = [];
            Data.(char(stns(i))).RFs_in_Time.ID.RFts(ind,:) = [];
            Data.(char(stns(i))).RFs_in_Time.ID.Errort(ind) = [];
            RF  = -Data.(char(stns(i))).RFs_in_Time.ID.RFs;
            TV =  Data.(char(stns(i))).RFs_in_Time.ID.Time;
            if (size(RF,1)) == 0;
                break
            end
                fig = figure
                plot(TV,RF)
                xlim([-5 20])
                title(char(stns(i)));
                datacursormode on  % Enable data cursor mode
                dcm_obj = datacursormode(fig);
                set(dcm_obj,'UpdateFcn',@myupdatefcn) % Set update function
                % Wait while the user to click
                disp('Click line to display a data tip, then press "Return"')
                pause
                info_struct = getCursorInfo(dcm_obj); % Export cursor to workspac
                close(fig)
        end
      end
    end
end
    
save(ofile,'Stns','Data');    
 