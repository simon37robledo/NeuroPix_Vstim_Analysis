j =1;
GoodRecordings =[40:43,49:54];%Anesthetized
Awake = [1:21, 28:36, 44:48];
amps = {};
sign = 0.005;
for ex =  GoodRecordings%GoodRecordingsPV%allGoodRec %GoodRecordings%GoodRecordingsPV%GoodRecordingsPV%selecN{1}(1,:) %1:size(data,1)
    %%%%%%%%%%%% Load data and data paremeters
    %1. Load NP class
    NP = loadNPclassFromTable(ex);

    p = NP.convertPhySorting2tIc(NP.recordingDir);
    label = string(p.label');
    goodU = p.ic(:,label == 'good');
    ampsAll = p.neuronAmp(label == 'good');
    respNeuronsMB = load(sprintf('pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponse;
    amps{j} = ampsAll(respNeuronsMB<=sign);

    j = j+1;

end
%%

Av = cell2mat(amps');
AvAwake = cell2mat(amps');

figure;histogram(Av);

% hold on
% histogram(AvAwake);
ylabel('# of neurons')
xlabel('Unit amplitudes (uV)')
title('Only neurons tuned to MB')
legend({'Anesthesia','Awake'})

