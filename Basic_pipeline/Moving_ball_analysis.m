path = '\\sil3\data\\Large_scale_mapping_NP\lizards\PV139\PV139_Experiment_6_2_24\Insertion1\catgt_PV139_Experiment_6_2_24_1_g0';

NP = NPAPRecording(path);

%%


%%%Moving ball inputs, NP, stimOn, stimOff, plots, neuron, ttlIndex


patternIndex = strfind(string(NP.recordingDir), "\catgt");

endIndex = patternIndex(1)-1;
stimDir = string(NP.recordingDir);
stimDir = extractBetween(stimDir,1,endIndex);

file = dir (stimDir);
filenames = {file.name};


file = dir (stimDir);
filenames = {file.name};
ballFiles = filenames(contains(filenames,"linearlyMovingBall"));
directions = [];
offsets = [];
sizes = [];
speeds = [];

j =1;
if size(ballFiles) ~= [0 0]

    for i = ballFiles
        ball= load(stimDir+"\"+string(i));

        directions = [directions cell2mat(ball.VSMetaData.allPropVal(17))];

        offsets = [offsets cell2mat(ball.VSMetaData.allPropVal(18))];

        sizes = [sizes cell2mat(ball.VSMetaData.allPropVal(19))];

        speeds = [speeds cell2mat(ball.VSMetaData.allPropVal(16))];

        direcNames = unique(directions);

        stimDurStats = cell2mat(ball.VSMetaData.allPropVal(38))*1000;
        interStimStats = cell2mat(ball.VSMetaData.allPropVal(28))*1000;

        j = j+1;
    end
    disp('Visual stats extracted!')
else
    disp('Directory does not exist!');
end


A = [stimOn' directions' offsets' sizes' speeds'];

C = sortrows(A,[2 3 4 5]);

%Sort directions:

directimesSorted = C(:,1)';

[Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round((stimDur+preBase)/bin));

%%%% Convolute matrix:

%%%1. Convolute in the 3rd dimension (trials)

[MrC]=ConvBurstMatrix(Mr,fspecial('gaussian',[1 5],3),'same');


[nT,nN,nB] = size(MrC);

trialDivision = nT/(length(unique(offsets))*length(unique(directions)));


%%% end of convolution



%%%%Create window of 2 offsets by 10 bins (500 ms) to scan
%%matrix and select highest mean druing scan.

% Define the size of the window
if a == 1 %PV67
    window_size = [4*trialDivision, 2]; %bin = 50 ms, 4 offsets
else
    window_size = [trialDivision, 2];
end


% Initialize the maximum mean value and its position %No
% need or concolutions.

max_position = zeros(nN,2);
max_mean_value = zeros(1,nN);

for u =1:nN
    % Slide the window over the matrix
    %unit matrix
    max_mean_value(u) = -Inf;
    for d = 1:length(unique(offsets))*trialDivision:nT
        uM = squeeze(MrC(d:d+length(unique(offsets))*trialDivision-1,u,:));
        for i = 1:size(uM, 1) - window_size(1) + 1
            for j = 1:size(uM, 2) - window_size(2) + 1
                % Extract the sub-matrix
                sub_matrix = uM(i:min(i+window_size(1)-1, end), j:min(j+window_size(2)-1,end));

                % Compute the mean value
                mean_value = mean(sub_matrix(:));

                % Update the maximum mean value and its position
                if mean_value > max_mean_value(u)
                    max_mean_value(u) = mean_value;
                    max_position(u,:) = [i+d-1, j];

                end


            end
        end
    end
end

spkRateR = max_mean_value;

epsilon = 0.01;

denom = mad(MbC,0)+epsilon; %mean(Mb,0)+epsilon; %

mSpk = mean(spkRateR);

Zscore = (spkRateR - (spkRateBM + spkRateR)/2)./denom;

RespU = GoodU_or(Zscore>2);


sDirections = sort(directions);

tunning = sDirections(max_position(:,1));



%MinusRB = (mean(spkRateR) -(mean(spkRateB)))./(mean(spkRateR) + (mean(spkRateB)));

MinusRB = spkRateR - spkRateBM;


%                 if a==2 & in ==1
%                     2+2
%
%                 end

if (a == 1 && in ==1) || (a == 2 && in ==1) || (a == 3 && in ==1)
    plotexamplesMB =0;

    if plotexamplesMB == 1
        %                     % Create direction arrows:
        %
        %                     % Define the x and y coordinates for the arrows
        %                     x = zeros(1,direcN);
        %                     y = categ/direcN:categ/direcN:categ;
        %                     theta = unique(directions);
        %
        %                     % Calculate the u and v components of the arrows
        %                     u = cos(theta);
        %                     v = sin(theta);
        %
        %                     % Normalize the arrow vectors to make them all the same length
        %                     arrow_length = 2;
        %                     u = u ./ sqrt(u.^2 + v.^2) * arrow_length;
        %                     v = v ./ sqrt(u.^2 + v.^2) * arrow_length;
        %
        %                     % Adjust the starting points of the arrows
        %                     x = x - u/2;
        %                     y = y - v/2;
        %
        %                     figd = figure;
        %                     h = quiver(x, y, u, v, 'AutoScale', 'off');
        %                     set(h, 'Color', 'k', 'LineWidth', 1);% Set the aspect ratio% Set the color to black and the line width to 2
        %                     set(gca, 'DataAspectRatio', [1 1 1]);
        %                     set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none'); % Remove ticks and axis borders
        %                     axis([-1 4 0 40]);
        %                     set(figd, 'Color', 'white');
        %                     print(figd, sprintf('%s-In.%d-MovBall-dir.png',expA{a},in),'-dpng');
        %
        %                     close all


        for u = eNeuron
            [Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round((stimDur+preBase*2)/bin));

            [nT,nN,nB] = size(Mr);

            fig = figure;
            categ = nT/trialDivision;
            offsetN = length(unique(offsets));
            direcN = length(unique(directions));
            imagesc(squeeze(Mr(:,u,:)));colormap(flipud(gray(64)));
            xline(preBase/bin,'k', LineWidth=1.5)
            xline(stimDur/bin+preBase/bin,'k',LineWidth=1.5)
            ylabel('Trials');xlabel('Time (ms)');
            m = max_mean_value(u) - spkRateBM(u);
            title(sprintf('U.%d-R.%.3f-B.%.3f-S.%.3f',u,max_mean_value(u),spkRateBM(u),Zscore(u)));
            xticks([0.5 (preBase/bin):10:nB])
            xticklabels([-preBase 0:10*bin:nB*bin])
            v = nT/direcN:nT/direcN:nT-1;
            yline(v+0.5,'r--', LineWidth=3);
            %                             hcb = colorbar();
            %                             title(hcb,'Spikes/sec');
            caxis([0 max(0.2,max(max_mean_value(u)))])
            hold on
            rectangle('Position', [max_position(u,2), max_position(u,1), window_size(2), window_size(1)], 'EdgeColor', 'r', 'LineWidth', 1.5,'LineStyle','-.');
            hold off
            prettify_plot
            %fig.Position = [-755   440   520   438];
            print(fig, sprintf('%s-In.%d-MovBall-dir-U%d-W%d-%dW.png',expA{a},in,u,window_size(1),window_size(2)),'-dpng');
            close

        end

    end
end



%end