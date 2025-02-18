

%%
entropies = zeros(1,length(respU));

for u = 1:length(respU)


    M = squeeze(RFuSTDir(:,:,:,u));

    % Find the maximum value and its index
    [maxValue, linearIndex] = max(M(:));

    % Convert the linear index to subscripts to find the slice
    [sliceIndex, ~, ~] = ind2sub(size(M), linearIndex);

    % Extract the corresponding 10x10 slice
    Mred = squeeze(M(sliceIndex, :, (size(M,3)-size(M, 2))/2+1: size(M,3)-(size(M,3)-size(M, 2))/2));

    %%Prepare matrix into the same style as rectGrid

    % Original 50x100 matrix
    %M = rand(50, 100); % Replace with your matrix

    % Number of blocks (offsetN)
    offsetN = 9;

    % Define edges for rows and columns
    rowEdges = round(linspace(1, size(Mred, 1) + 1, offsetN + 1));
    colEdges = round(linspace(1, size(Mred, 2) + 1, offsetN + 1));

    % Initialize the resulting matrix
    reducedMatrix = zeros(offsetN, offsetN);

    % Compute mean for each block, omitting NaNs
    for i = 1:offsetN
        for j = 1:offsetN
            % Extract block using the calculated edges
            block = Mred(rowEdges(i):rowEdges(i+1)-1, colEdges(j):colEdges(j+1)-1);
            % Compute mean, omitting NaNs
            reducedMatrix(i, j) = mean(block(:), 'omitnan');
        end
    end

    % Normalize to create a probability distribution
    reducedMatrix = reducedMatrix / sum(reducedMatrix(:));

    % Convert to an image-like format and calculate entropy
    % (scale to [0, 1] for compatibility with `entropy`)
    P_scaled = mat2gray(reducedMatrix);
    figure;imagesc((P_scaled));
    set(gca, 'YDir', 'reverse', 'XDir', 'reverse');

    entropies(u) = entropy(P_scaled);
end


cd(NP.recordingDir)
sign = '0.005';

save(sprintf('Entropies-MB-RF-respU-%s-%s.mat',sign,NP.recordingName),'entropies')


