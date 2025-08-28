function[RF_Time,optimal_alpha] = ETMTM_new(P,D,time,TB,NT,t0,t1,Faza)
        %TB = 4; NT = 7;
        % Findings: when P arrival is not in the center of the window, the
        % amplitudes are not unity at the beginning and decreasing from there on.
        % Instead they peak at the time shift which corresponds to the middle index
        % in the P time window.
        %P = P (vertical)
        %D = SV (raidal)
        % TB  = time bandwidth product (usually between 2 and 4)
        % NT  = number of tapers to use, has to be <= 2*TB-1
        %
        % TB = 4; NT = 7; %choice of TB = 4, NT = 3 is supposed to be optimal
        % t0 = -5; t1 = max(time); % This defines the beginning and end of the lag
        % times for the receiver function
        
        % Ved wrote MTM for MATLAB, which has the added advantage of
        % finding the optimal damping parameter.
        
        % Find sampling interval
        dt = time(2) - time(1);
        
        % Flip time axis in case of Sp
        if(strcmp(Faza,'Sp'))
            D = fliplr(D); P = fliplr(P);
%             tmp0 = t0;    tmp1 = t1;
%             t0   = -tmp1; t1   = -tmp0;
        end
        
        % Length of moving time window in seconds
        % THE ARRIVAL TIME OF THE PARENT PHASE HAS TO BE IN THE MIDDLE OF THE FIRST
        % WINDOW. ALSO TRUE FOR SP WHERE THE WAVEFORMS ARE FLIPPED
        
        win_len = (length(D)*dt)/2; % this window length is computed with 
                                    % the assumption that Auto_Prep prepared
                                    % the waveforms to have 1/4 of the time
                                    % before the arrival and 3/4 after (or 
                                    % vice versa for Sp). For signal that
                                    % starts at 1/4 to be in the middle of
                                    % the first window, it has to be length
                                    % 1/2.
        Nwin = round(win_len/dt);
        
        % Fraction of overlap overlap between moving time windows. As your TB
        % increases, the frequency smearing gets worse, which means that the RFs
        % degrate at shorter and shorter lag times. Therefore, as you increase TB,
        % you should also increase Poverlap.
        Poverlap = 0.99;
        
        % Create moving time windowed slepians
        starts = 1:round((1-Poverlap)*Nwin):length(P)-Nwin+1;
        
        
        % Length of waveforms;
        nh = length(P);
        
        % Construct Slepians
        [Etmp,lambdas] = dpss(Nwin,TB);
        
        E = zeros(length(starts)*size(Etmp,2),nh);
        n = 0;
        
        NUM = zeros(size(P));  DEN = zeros(size(D));
        
        for j = 1:length(starts)
            for k = 1:NT
                n = n + 1;
                E(n,starts(j):starts(j)+Nwin-1) = transpose(Etmp(:,k));
                
                tmp1 = fft(E(k,:).*P); % always stick to first time window
                tmp2 = fft(E(n,:).*D); % allow moving time windows
                
                NUM = NUM + lambdas(k)*conj(tmp1).*tmp2;
                DEN = DEN + lambdas(k)*conj(tmp1).*tmp1;
                
            end
        end
        
        % Search range for optimal damping parameter alpha
        alphas = logspace(-6,6,65)*var(D(round(end/4):3*round(end/4)))*length(P);
        
        % WE MULTIPLY ALPHA BY THE VARIANCE OF THE MIDDLE HALF OF THE DAUGHTER
        % RECORD AND BY THE LENGTH OF THE WAVEFORM IN POINTS (DUE TO MATLAB
        % NORMALIZATION).
        
        % Initialize misfit, RF magnitude, and RF structures
        misfit = zeros(size(alphas));
        magntd = zeros(size(alphas));
        
        % Now, calculate misfit and RF size for each alpha
        for kj = 1:length(alphas)
            
            % Calculate  RF
            tmp = real(ifft(NUM./(DEN + alphas(kj))));
            
            % Normalize  RF
            nrm = max(real(ifft(NUM./(NUM + alphas(kj)))));
            RF_Time = reshape(fftshift(tmp./nrm),[1 nh]);
            
            % Calculate the misfit between predicted and observed daughters
            tmp = ifft(conj(fft(RF_Time)).*fft(P));
            misfit(kj) = nansum(abs(D - tmp));
            magntd(kj) = nansum(abs(tmp));
        end
        
        % Find optimal alpha - use a 7 point mask to smooth and find
        % the gradient all in one step in order to supress the 
        % meaningless oscillations in the misfit due to noise 
        tmp = imfilter(log((misfit./std(misfit)).^2+(magntd./std(magntd)).^2),...
                       -transpose([gradient(gausswin(7))]),'replicate'); 
        %plot(log(alphas),tmp,'rs'); hold on; 
        %plot(log(alphas),log((misfit./std(misfit)).^2+(magntd./std(magntd)).^2));
        % Then, select the alpha at which the misfit grows most quickly. 
        [~,j2] = max(tmp); optimal_alpha = alphas(j2); 

        % Calculate optimal RF
        tmp = real(ifft(NUM./(DEN + optimal_alpha)));
        
        % Normalize optimal RF
        nrm = max(real(ifft(NUM./(NUM + optimal_alpha))));
        RF_Time = reshape(fftshift(tmp./nrm),[1 nh]);
          
        % If Sp, then flip back the receiver function
        
        if(strcmp(Faza,'Sp')), RF_Time = fliplr(RF_Time); end;
           % Interpolate to desired
        RF_Time = interp1(dt*(-nh/2+0.5:1:nh/2-0.5),RF_Time,t0:dt:t1,'linear',0);
      
        
     %  display(['Optimal value of damping = ' num2str(alphas(j2)) '...']);
        
    end