% Example data: Raster matrix with 200 trials and 1000 time bins (1 ms per bin)

u = 57;
cd(NP.recordingDir)

%         mr = squeeze(BuildBurstMatrix(goodU(:,u),round(p.t),round(directimesSorted),round(stimDur)));

        trialsPerAngle = trialDivision*offsetN*speedN*sizeN*orientN;
        
        ZscoreRaster = load(sprintf('ZscoreRaster-%d-%s',N_bootstrap,NP.recordingName)).ZscoreRaster;
        ZscoreRaster = squeeze(ZscoreRaster(:,u,:));

        [nT,nB] = size(ZscoreRaster);
% 
%         mr = max(squeeze(ZscoreRaster(:,u,:)),[],2); %%Take mean per offset, then take max offset as result to plot
%         
        trialsPerOffset = trialDivision*sizeN; 

        maxZ = zeros(1,nT/trialsPerOffset);
        maxZpos = zeros(1,nT/trialsPerOffset);
        j=1;
       
        for i =1:trialsPerOffset:nT
            meanZSperOffset = median(ZscoreRaster(i:i+trialsPerOffset-1,:),1);

            [maxZ(j) maxZpos(j)] = max(meanZSperOffset);

            j =j+1;
        end

      % 

        % Initialize arrays to store spike rates and SEM for each angle
        tuningValZS = zeros(1, direcN); % 4 angles (0, 90, 180, 270)
        sem_values = zeros(1, direcN);   % 4 SEM values
        duration = 300;

        % Loop through each angle and calculate the max ZS and SEM
        for i = 1:direcN
            trials_for_angle = find(uDir(i) == C(:,2)); % Find trials for each angle (0, 90, 180, 270)

            %Calculate max z-score per offsets within a direction 
            [tuningValZS(i) trPosition]= max(maxZ(i*offsetN-(offsetN-1):max(trials_for_angle)/trialDivision));

            %%Obtain trialDivision (e.g. 10) trials from bin position of max
            %%Z-score within direction
            TrialValsOfMaxZS = ZscoreRaster(trPosition*trialDivision-(trialDivision-1):trPosition*trialDivision,maxZpos(trPosition));
            
%             % Calculate 
%             trial_spike_rates = max(mr(trials_for_angle, :), 2) / (stimDur/1000); % Spike rate per trial (spikes per second)
% 
            % Calculate the mean spike rate and SEM
            sem_values(i) = std(TrialValsOfMaxZS) / sqrt(length(TrialValsOfMaxZS));
        end


        % Plot the tuning curve with error bars
        fig = figure;
        angles_deg = rad2deg(uDir); % Angles in degrees
        bar(angles_deg, tuningValZS, 'FaceColor', 'k','FaceAlpha',0.5); % Bar plot for spike rates
        hold on;
        errorbar(angles_deg, tuningValZS, sem_values, 'k', 'LineStyle', 'none', 'LineWidth', 1.5); % Error bars
        xlabel('Stimulus Angle (degrees)');
        ylabel('Z-score');