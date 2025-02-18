baseline = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-interStimStats)/bin),round((interStimStats)/bin));
baseline = single(baseline);
[nT,nN,nB] = size(baseline);
% Bootstrapping settings
N_bootstrap = 1000; % Number of bootstrap iterations
boot_means = zeros(N_bootstrap, nN,'single');
resampled_indicesTr = single(randi(nT, [nT, N_bootstrap]));% To store bootstrapped means
resampled_indicesTi = single(randi(nB, [nB, N_bootstrap]));

kernel = ones(trialDivision, duration) / (trialDivision * duration); % Normalize for mean
% Start a parallel pool (if not already started)
if isempty(gcp('nocreate'))
    parpool; % Start a pool with the default number of workers
end

tic
parfor i = 1:N_bootstrap
    % Resample trials with replacement
    resampled_trials = baseline(resampled_indicesTr(:, i), :,resampled_indicesTi(:, i));
    for ui = 1:nN
        % Extract the slice for the current unit (t x b matrix)
        slice = resampled_trials(:, ui, :);
        slice = squeeze(slice); % Result is t x b

        % Compute the mean using 2D convolution
        means = conv2(slice, kernel, 'valid'); % 'valid' ensures the window fits within bounds

        % Find the maximum mean in this slice
        boot_means(i, ui) = max(means(:));
    end
end
toc
%%


%%%Generate indexes so that the convolutions doesn't mix between slices of
%%%neurons.
b = nB:nB:nB *nN;
% Generate the relative indices for each range
relative_indices = (0:(nB - duration)) - (nB - 1); % Creates [-nB+1, ..., -duration+1]

% Use broadcasting to compute all ranges
ranges = b.' + relative_indices; % A matrix where each row is a range
vec = ranges(:).';

baselineResh = permute(baseline, [1, 3, 2]);

boot_meansConv = zeros(N_bootstrap, nN,'single');

% Start a parallel pool (if not already started)
if isempty(gcp('nocreate'))
    parpool; % Start a pool with the default number of workers
end

tic

parfor i = 1:N_bootstrap
    resampled_trials = baselineResh(resampled_indicesTr(:, i), resampled_indicesTi(:, i), :);

    A_reshaped = reshape(resampled_trials, nT, nB * nN); % Collapse u into the second dimension

    % Apply 2D convolution across t and b
    conv_result = conv2(A_reshaped, kernel, 'valid'); % Result is (t - window_t + 1) x (b*u - window_b + 1)
    
    conv_resultTr = conv_result(:,vec);

    % Reshape the result back to separate u
    conv_result_reshaped = reshape(conv_resultTr, nT-trialDivision+1, nB - duration + 1, nN); % Back to (t - window_t + 1) x (b-window_b+1) x u

    % Compute the maximum mean for each u
    max_means = squeeze(max(conv_result_reshaped, [], [1, 2])); % Max over t and b for each u
    
    boot_meansConv(i,:) = max_means;
end
toc
%%
boot_meansW = zeros(N_bootstrap, nN,'single');
% Initialize variable to store the maximum mean
max_mean = -inf; % Start with a very small value
% Perform bootstrapping
for i = 1:N_bootstrap
    % Resample trials with replacement
    resampled_trials = baseline(resampled_indicesTr(:, i), :,resampled_indicesTi(:, i));
    for w = 1:(nT - trialDivision + 1)  % Loop over rows
        for j = 1:(nB - duration + 1)  % Loop over columns
            % Extract the current window
            window = resampled_trials(w:w+trialDivision-1, j:j+duration-1);

            % Compute the mean of the current window
            window_mean = mean(window(:));

            % Update the maximum mean if the current window mean is larger
            max_mean = max(max_mean, window_mean);
        end
    end
    % Compute mean firing rate across resampled trials
    boot_meansW(i, :) = max_mean;%mean(resampled_trials, [1 3]);
end
toc
bootSections = 50; %%%Off course verything is going to be signifficant because the eman is in different dimensions
%%

% Compute confidence intervals
lower_bound = prctile(boot_means, 2.5, 1); % 2.5th percentile (lower CI)
upper_bound = prctile(boot_means, 97.5, 1); % 97.5th percentile (upper CI)

[respVal respVali]= max(NeuronVals(:,:,1),[],2);

 Mr = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted)/bin),round((stimDur+duration)/bin)); %response matrix

%%% Calculate p-value & Filter out neurons in which max response window is empty for more than
%%% 60% of trials

pvalsResponse = zeros(1,nN);

for u = 1:nN
    posTr = NeuronVals(u,respVali(u),2);
    posBin = NeuronVals(u,respVali(u),3);

    maxWindow = squeeze(Mr(posTr*trialDivision-trialDivision+1:posTr*trialDivision,u,posBin:posBin+duration-1));

    emptyRows = sum(all(maxWindow == 0, 2));

    pvalsResponse(u) = mean(boot_means(:,u)>respVal(u));

    if emptyRows/trialDivision > 0.6
        pvalsResponse(u) = 1;
    end

end

figure;imagesc(maxWindow)
figure;imagesc(squeeze(Mr(:,u,:)));
rectangle('Position', [posBin,posTr*trialDivision-trialDivision+1, duration, trialDivision],...
                        'EdgeColor', 'r', 'LineWidth', 1.5,'LineStyle','-.');


%%
% Example for a specific neuron and time bin
neuron_idx = 1; % Index of neuron
bin_idx = 10;   % Index of time bin
disp(['95% CI for neuron ' num2str(neuron_idx) ', bin ' num2str(bin_idx) ':']);
disp(['Lower Bound: ' num2str(lower_bound(1, neuron_idx, bin_idx))]);
disp(['Upper Bound: ' num2str(upper_bound(1, neuron_idx, bin_idx))]);

% Visualization for a specific neuron and time bin
figure;
histogram(squeeze(boot_means(:, neuron_idx, bin_idx)), 30); % Distribution of bootstrapped means
xlabel('Mean Firing Rate');
ylabel('Frequency');
title(['Bootstrap Distribution (Neuron ' num2str(neuron_idx) ', Bin ' num2str(bin_idx) ')']);
hold on;
xline(lower_bound(1, neuron_idx, bin_idx), 'r', 'Lower CI');
xline(upper_bound(1, neuron_idx, bin_idx), 'r', 'Upper CI');
hold off;
