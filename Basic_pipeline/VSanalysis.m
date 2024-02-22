%%%
path = '\\sil3\data\\Large_scale_mapping_NP\lizards\PV139\PV139_Experiment_6_2_24\Insertion1\catgt_PV139_Experiment_6_2_24_1_g0';

NP = NPAPRecording(path);



%%

% function []=VSanalysis(NP,stimDir,varargin)

stimDir = '\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\PV67\PV67_experiment_5_7_23\Insertion1';

%obtain stim file names
file = dir (stimDir);
filenames = {file.name};
stimFiles = filenames(contains(filenames,"linearlyMovingBall"|"fullFieldFlash"|"rectGrid"|"StaticDriftingGrating"));

%Extract time of stimulus presentation from files
for i = 1:length(stimFiles)
    numbers = regexp(stimFiles{i},'\d+','match');
    concatNum = strcat(numbers{:});
    orderStim(i) = str2double(concatNum);     
end

%Sort stim files acordying to the time they were presented
[s,i] = sort(orderStim);
ttlindex = [1:length(stimFiles);i];

%Order stim files
for i = 1:length(stimFiles)
    stimFilesOrdered{i} = stimFiles{ttlindex(2,i)};
end

%Operate stim functions

if allStim ==1
    for i = 1:stimFiles

        stimStats = load(stimFilesOrdered{i}); 


       

    end

else

end












