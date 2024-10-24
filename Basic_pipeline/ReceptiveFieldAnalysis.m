%% RF plots mov ball/rectGrid

%% Do the sumary graphs. Boxplot for moving ball. 

cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

examplesSDG =[1 2 3 4 5 6 7 8 9 10 11 12 13 14 29 30 31 32 40 41 42 43];

%%
for ex = 41

     path = convertStringsToChars(string(data.Base_path(ex))+filesep+string(data.Exp_name(ex))+filesep+"Insertion"+string(data.Insertion(ex))...
        +filesep+"catgt_"+string(data.Exp_name(ex))+"_"+string(data.Insertion(ex))+"_g0");
    try %%In case it is not run in Vstim computer, which has drives mapped differently
        cd(path)
    catch
        originP = cell2mat(extractBetween(path,"\\","\Large_scale"));
        if strcmp(originP,'sil3\data')
            path = replaceBetween(path,"","\Large_scale","W:");
        else
            path = replaceBetween(path,"","\Large_scale","Y:");
        end
        cd(path)
    end
    NP = NPAPRecording(path);

     %2. Extract moving ball statistics
    patternIndex = strfind(string(NP.recordingDir), "\catgt");

    endIndex = patternIndex(1)-1;
    stimDir = string(NP.recordingDir);
    stimDir = extractBetween(stimDir,1,endIndex);

    file = dir (stimDir);
    filenames = {file.name};

    file = dir (stimDir);
    filenames = {file.name};
    ballFiles = filenames(contains(filenames,"linearlyMovingBall"));

    
        if isempty(ballFiles)
            %disp()
            w= sprintf('No moving ball files where found in %s. Skipping into next experiment.',NP.recordingName);
            warning(w)
            continue
        end

    directions = [];
    offsets = [];
    sizes = [];
    speeds = [];
    orientations = [];

    j =1;

    if size(ballFiles) ~= [0 0]

        for i = ballFiles %Extract stim properties
            ball= load(stimDir+"\"+string(i));

            if ~isempty(find(strcmp(ball.VSMetaData.allPropName,'orientations'))) %Check if orientations are present (grid inside moving object).
                orientations = [orientations cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'orientations'))))];
            else %No orientations property exist
                orientations = [orientations zeros(1,cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'nTotTrials')))))];
            end
            if isempty(orientations) %orientations property exists but is empty
                orientations = [orientations zeros(1,cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'nTotTrials')))))];
            end

            directions = [directions cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'directions'))))];

            offsets = [offsets cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'offsets'))))];

            sizes = [sizes cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'ballSizes'))))];

            speeds = [speeds cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'speeds'))))];

            interStimStats = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'interTrialDelay'))))*1000;

            j = j+1;
        end
        disp('Visual stats extracted!')
    else
        disp('Directory does not exist!');
    end

    if ~isempty(find(strcmp(ball.VSMetaData.allPropName,'orientationFreq')))
        Freq = [ cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'orientationFreq'))))];
    end

    
    uDir = unique(directions);
    uSize = unique(sizes);
    offsetN = length(unique(offsets));
    direcN = length(unique(directions));
    speedN = length(unique(speeds));
    sizeN = length(unique(sizes));
    orientN = length(unique(orientations));
    nT = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'nTotTrials'))));
    trialDivision = nT/(offsetN*direcN*speedN*sizeN*orientN); %Number of trials per unique conditions
    categ = nT/trialDivision;



RFu = load(string(path)+filesep+"RFuC-"+string(NP.recordingName)).RFu;

pvalTi= load(sprintf('pvalTime-%s',NP.recordingName)).pvalTi;

respU = find(pvalTi <0.02);


% %%%%%%% Find max receptive field and how localized it is (entropy): It runs
%%%%%%% a ball shaped window that moves across the image of the previously
%%%%%%% calculated convolution. 
RFuSizeSav = load(sprintf('RFuSizeC-%s',NP.recordingName)).RFuSize;
RFuDirSav = load(sprintf('RFuDirC-%s',NP.recordingName)).RFuDir; 
NormVideo = load(sprintf('NormVideo-%s',NP.recordingName)).NormVideo; 

sizeN = size(RFuSizeSav,1);
xSize= size(RFuSizeSav,2);
ySize= size(RFuSizeSav,1);

%Size analysis -> 1. hacer la grafica que pide Mark, diferentes Neuronas.
%2. Calcular entropia en experimento y graficar contra la profundidad
coorRect = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'rect'))));

convMeanSize = zeros(sizeN,coorRect(4)/10,coorRect(3)/10,length(respU));
maxValdelayS = zeros(sizeN,length(respU));
maxInddelayS = zeros(sizeN,length(respU));



sumVid = squeeze(sum(NormVideo,3));
repeatedNormvid = single(repmat(reshape(sumVid,size(sumVid,1),size(sumVid,2),1,size(sumVid,3)),[1 1 size(RFu,3) 1]));
%repeatedNormvid(repeatedNormvid==0) = 1;

% implay(A(:,:,:,9))
%
NormVideo2= NormVideo;
NormVideo2(NormVideo2==0) =1;
%%
%figure;imagesc(squeeze(repeatedNormvid(:,:,5,9)));colorbar %do norm of baseline at the end, %change normalization, %check eye videos of pv35
for s=1:sizeN

    % Parameters for the circular mask
    radius = uSize(s)/10; % Define the radius of the circular window
    % Create a circular mask
    [X, Y] = meshgrid(-radius:radius, -radius:radius);
    circular_mask = (X.^2 + Y.^2) <= radius^2;

    % Normalize the circular mask so that the sum of the mask elements is 1
    circular_mask = single(circular_mask);
    circular_mask = circular_mask / sum(circular_mask(:));
    
    A = double(squeeze(RFuSizeSav(s,:,:,:,:)))./double(NormVideo);

    %  Old method:   % Use convolution to find the sum of elements inside the circular window
    %     A = squeeze(RFuSizeSav(s,:,:,:,:))./NormVideo;
    %     A_reshaped = reshape(A, size(A, 1), size(A, 2), []);
    %
    %     convolved_result = convn(A_reshaped, circular_mask, 'same');
    %     % Reshape the convolved result back to the original 4D structure
    %     convolved_result_reshaped = reshape(convolved_result, size(A, 1), size(A, 2), size(A, 3), size(A, 4));
    %
    %     padding =150;
    %     A_padded = padarray(A, [0, 0, padding,0], 'post');
    %
    %     % Note: The filter should be expanded to 3D by adding a singleton dimension
    %     circular_mask_3d = reshape(circular_mask, [size(circular_mask, 1), size(circular_mask, 2), 1, 1]);
    %
    %     % Perform the convolution
    %     convolved_result_reshaped = convn(A_padded, circular_mask, 'same');
    %
    %
    %     % Initialize arrays to store results
    %     max_means = zeros(size(convolved_result, 3), size(convolved_result, 4));
    %     center_row = zeros(size(convolved_result, 3), size(convolved_result, 4));
    %     center_col = zeros(size(convolved_result, 3), size(convolved_result, 4));
    %
    %     % Find the maximum values and their positions for each slice across the 3rd and 4th dimensions
    %     for j = 1:size(convolved_result_reshaped, 4)
    %         for i = 1:size(convolved_result_reshaped, 3)
    %             % Extract the 2D slice
    %             slice = convolved_result_reshaped(:, :, i, j);
    %             % Find the maximum value and its position
    %             [max_means(i, j), max_idx] = max(slice(:));
    %             [max_row, max_col] = ind2sub(size(slice), max_idx);
    %             % Adjust positions to get the center of the circular mask
    %             center_row(i, j) = max_row + radius;
    %             center_col(i, j) = max_col + radius;
    %         end
    %     end
    %
    %     [max_across_3rd, idx_3rd] = max(max_means, [], 1);
    %
    %     maxValdelayS(s,:) = max_across_3rd;
    %
    %     maxInddelayS(s,:) = idx_3rd;


    %     for u = 1:length(respU)
    %         convMeanSize(s,:,:,u) =  squeeze(convolved_result_reshaped(:,:,idx_3rd(u),u));
    %     end
    for t = 1:size(A, 4) %For neuron
        entropyValues = zeros(1, size(A, 3));  % To store entropy values for each slice in 3rd dimension

        % Preallocate matrix to store circular mean results
        meanSlices = zeros(size(A, 1), size(A, 2), size(A, 3));

        % Loop over the 3rd dimension (depth, z slices)
        for k = 1:size(A, 3) %for frame
            % Perform convolution in the first two dimensions with the circular mask
            meanSlices(:,:,k) = conv2(A(:,:,k,t), circular_mask, 'same');

            % Compute the entropy of the resulting 2D slice
            currentSlice = meanSlices(:,:,k);
            entropyValues(k) = entropy(currentSlice);  % MATLAB's entropy function
        end

        % Find the index of the slice with the lowest entropy in the 3rd dimension
        [~, minEntropyIdx] = min(entropyValues(round(size(A,3)/2):end));

        % Store the slice with the lowest entropy for this 4th dimension (time, etc.)
        convMeanSize(s,:,:,t) = meanSlices(:,:,minEntropyIdx+round(size(A,3)/2)-1);
        
        maxInddelayS(s,t) = minEntropyIdx+round(size(A,3)/2)-1;

    end %end of looping neurons
end %End of loopinf size

size(RFu)

% sumVid = squeeze(sum(NormVideo,3));
% 
% figure;imagesc(squeeze(sumVid(:,:,9)))

A_normalized = squeeze(convMeanSize(1,:,:,:));

% implay(squeeze(A(:,:,:,9)))
% implay(squeeze(NormVideo(:,:,:,9)))
% 
max(A_normalized(:,:,8),[],'all')
min(A_normalized(:,:,8),[],'all')

figure;imagesc(flip(squeeze(A_normalized(:,:,8))));
figure;imagesc(flip(squeeze(A(:,:,90,9))));
figure;imagesc(flip(squeeze(meanSlices(:,:,90))));

implay(NormVideo)


for k = 1:length(respU)
    slice = A_normalized(:, :, k);
    min_val = min(slice(:));
    max_val = max(slice(:));
    A_normalized(:, :, k) = (slice - min_val) / (max_val - min_val);
end

[m i] = max(max_means(:,27))
figure;plot(max_means(:,27)');xline(90);xline(i,'Color','r');xlabel('frames');ylabel('max mean value')
implay(A_normalized)

end





