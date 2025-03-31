%%%%%%% tuning analysis




cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

GoodRecordingsPV =[40:43 49:54];

recordings = GoodRecordingsPV;

takeMedian = 1;

tuningCurve = cell(1,numel(recordings));
SEM = cell(1,numel(recordings));
DSI = cell(1,numel(recordings));
OSI =cell(1,numel(recordings));
Theta = cell(1,numel(recordings));
preferDir = cell(1,numel(recordings));

MBc =[];
i =1;
animal =0;
for ex =  GoodRecordingsPV
    %%%%%%%%%%%% Load data and data paremeters

    %1. Load NP class
    path = convertStringsToChars(string(data.Base_path(ex))+filesep+string(data.Exp_name(ex))+filesep+"Insertion"+string(data.Insertion(ex))...
        +filesep+"catgt_"+string(data.Exp_name(ex))+"_"+string(data.Insertion(ex))+"_g0");
    try %%In case it is not run in Vstim computer, which has drives mapped differently
        cd(pathE)
    catch
        try
            originP = cell2mat(extractBetween(path,"\\","\Large_scale"));
            if strcmp(originP,'sil3\data')
                path = replaceBetween(path,"","\Large_scale","W:");
            else
                path = replaceBetween(path,"","\Large_scale","Y:");
            end

            cd(path)
        catch

            if strcmp(originP,'sil3\data')
                path = replaceBetween(path,"","\Large_scale","\\sil3\data");
            else
                path = replaceBetween(path,"","\Large_scale","\\sil1\data");
            end
            cd(path)

        end
    end
    NP = NPAPRecording(path);

    N_bootstrap =1000;
    cd(NP.recordingDir)
    respNeuronsMB = load(sprintf('pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponse;
    sign = 0.05; %%%Significance level used to calculate receptive fields
    respU = find(respNeuronsMB<sign);
        
    tempTC = load('tuningCurveAllOffsets').tuningCurve;;

    if size(tempTC,2) > 4
        tempTC = tempTC(:,1:2:size(tempTC,2)); %%Select only up, left, down, and right (0,90,180,270) if the recording has more than 4 directions
    end


    tuningCurve{i} = tempTC;


    DSIt= load(sprintf('Direction-Selectivity-Index-%s',NP.recordingName)).DSI;
    DSIt = DSIt(respNeuronsMB<sign);
    DSI{i} = DSIt;
    
    OSIt = load(sprintf('Orientation-Tuning-Index-%s',NP.recordingName)).L;
    OSIt = OSIt(respNeuronsMB<sign);
    OSI{i} = OSIt;
    
    Thetat = load(sprintf('Angle-prefer-%s',NP.recordingName)).pI;
    Thetat = Thetat(respNeuronsMB<sign);
    Theta{i} = Thetat;
    [maxVal preferD] = max(tempTC(respNeuronsMB<sign,:),[],2);

    preferDir{i} = preferD;

     %%%Add color vectors accordying to animals
    if string(data.Animal_ID{ex}) ~= string(data.Animal_ID{ex-1}) %wont work if you start with the first animal (noisy animal)
        animal = animal+1;
        animalName{animal} = data.Animal_ID{ex};
    end

    MBc = [MBc ...
            zeros(1,length(DSI{i}))+ animal]; %%%Add animal number for colors

    i=i+1;

end

%% plot histogram with a histogram per direction that shows the distribution of OSI and DSI values

cats = categorical([ones(1,length(cell2mat(DSI'))) ones(1,length(cell2mat(DSI')))+1 ones(1,length(cell2mat(DSI')))+2])';

fig = figure;
tiledlayout(2,1, "TileSpacing","loose","Padding","loose")

%%%%%%%%%%%
nexttile %%% all

swarmchart(cats, [cell2mat(DSI');cell2mat(OSI');cell2mat(DSI')-cell2mat(OSI')], 10, 'filled','MarkerFaceAlpha',0.5)

set(gcf,'Color','w');%
xticklabels({'DSI','OSI','DSI-OSI'})
ax = gca; % Get current axis
ax.XAxis.FontSize = 7; % Set font size of x-axis tick labels
ylabel('Tuning strength')
yline(0,'LineWidth',1,'Color','k')

% %%%%%%%%%%
% nexttile %%% dsi > 0.7
% 
% thres1 = 0.7;
% vDSI = cell2mat(DSI');
% 
% vOSI = cell2mat(OSI');
% 
% 
% diffDSI_OSI = vDSI(vDSI>thres1)-vOSI(vDSI>thres1);
% 
% val = [vDSI(vDSI>thres1);vOSI(vDSI>thres1);diffDSI_OSI];
% 
% cats = categorical([ones(1,sum(vDSI>thres1)) ones(1,sum(vDSI>thres1))+1 ones(1,sum(vDSI>thres1))+2])';
% 
% colors = [MBc(vDSI>thres1)'; MBc(vDSI>thres1)';MBc(vDSI>thres1)'];
% 
% swarmchart(cats, val, 10,'filled','MarkerFaceAlpha',0.8)
% 
% yline(thres1,'LineWidth',2,'Label',string(0.7),'LabelHorizontalAlignment','right','FontSize',7)
% xticklabels({'DSI > 0.7','OSI','DSI-OSI'})
% ax = gca; % Get current axis
% ax.XAxis.FontSize = 7; % Set font size of x-axis tick labels
% ylabel('Tuning strength')


% %%%%%%%
% nexttile %%% osi > 0.7
% 
% 
% thres1 = 0.7;
% 
% 
% diffOSI_DSI = vOSI(vOSI>thres1)-vDSI(vOSI>thres1);
% 
% val = [vDSI(vOSI>thres1);vOSI(vOSI>thres1);diffOSI_DSI];
% 
% cats = categorical([ones(1,sum(vOSI>thres1)) ones(1,sum(vOSI>thres1))+1 ones(1,sum(vOSI>thres1))+2])';
% 
% colors = [MBc(vOSI>thres1)'; MBc(vOSI>thres1)';MBc(vOSI>thres1)'];
% 
% swarmchart(cats, val, 10,'filled','MarkerFaceAlpha',0.8)
% 
% yline(thres1,'LineWidth',2,'Label',string(0.7),'LabelHorizontalAlignment','right','FontSize',7)

% set(gcf,'Color','w');%
% % xticklabels({'DSI','OSI > 0.7','OSI-DSI'})
% ax = gca; % Get current axis
% ax.XAxis.FontSize = 7; % Set font size of x-axis tick labels

%fig.Position = [1269         362         189         305];

ylabel('Tuning strength')
ylim([-0.1 1])
cd('\\sil3\data\Large_scale_mapping_NP\Figs paper\1stFigure')

nexttile
plot(vDSI,vOSI,'.');
xlabel('DSI')
ylabel('OSI')
axis equal
xlim([0 1])
ylim([0 1])
hold on
plot([0,1],[0,1],'LineWidth',1,'Color','k')
title('Median per direction')
print(fig, 'tuningIndexesMovBall.pdf', '-dpdf', '-r300', '-vector');

%% Histogram of prefered angles

%%Create a round histogram divided into 4 direction. Within each direction,
%%DSI values are sorted

angles = cell2mat(preferDir');
DSIv = cell2mat(DSI');
OSIv = cell2mat(OSI');

figure;swarmchart(categorical(angles),OSIv,10, 'filled','MarkerFaceAlpha',0.5);

set(gcf,'Color','w');%
ax = gca; % Get current axis
ax.XAxis.FontSize = 7; % Set font size of x-axis tick labels
ylabel('DSI')
yline(0,'LineWidth',1,'Color','k')

fig.Position = [1269         362         189         305];


fig = figure; tiledlayout(2,2)









