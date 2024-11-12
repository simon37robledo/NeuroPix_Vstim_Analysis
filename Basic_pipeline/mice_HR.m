%Mice HR analysis

%Load recording and create neuropixel class

base_dir = "\\sil3\data\Mice_exp\mouse1\Mice_exp_29_11_23\Insertion1\catgt_Mice_exp_28_11_23_1_g0";
%base_dir = "\\sil3\data\Mice_exp\mouse1\Mice_exp_29_11_23\Insertion1\catgt_Mice_exp_29_11_23_1_g0";%change to NP.recordingDir
run_dir = "\\132.66.45.127\data\Large_scale_mapping_NP\SpikeGLX\"; %Folder that has catgt and tprime subfolders
insertion = "1";
t1 = 0.8;
t2 = 0.9;
dur = 12;
dur_var = 5;
ttl_index =-1;
%%
recPath1 = '\\132.66.45.127\data\Large_scale_mapping_NP\Mice_experiments\Mouse_rightV1_NP_8_7_24\Insertion1\catgt_Mouse_rightV1_NP_8_7_24_1_g0'; 
%recPath1 = '\\132.66.45.127\data\Large_scale_mapping_NP\Mice_experiments\mouse1\Mice_exp_28_11_23\Insertion1\catgt_Mice_exp_28_11_23_1_g0'; 
NP1= NPAPRecording(recPath1);
%%
cd(NP1.recordingDir)

niStartsample = 50867337;
apStartsample = 144056030;

calSFheadstage = 30000.37799442897;
calSFnidaq = 10593.236768802228;

niStartTime = niStartsample/calSFnidaq;
apStartTime = apStartsample/calSFheadstage;
startDiff = niStartTime-apStartTime;
totalSamplesNi = round(3365.0159792817526*calSFnidaq);
totalSamplesAP = round(3365.2506*calSFheadstage); %Modify later with accurate probe info.

samplesTimeAp = (1:totalSamplesAP)/calSFheadstage;
samplesTimeNi = (niStartsample:totalSamplesNi)/calSFnidaq;

samplesTimeNi(end)-startDiff;
samplesTimeAp(end)



figure;
plot((niStartsample:totalSamplesNi),samplesTimeNi); hold on; plot((apStartsample:totalSamplesAP),samplesTimeAp); 

figure;
plot(samplesTimeAp-samplesTimeNi)


%[astim, astimon, astimoff, astimdur,astiminter] = time_diode_synced(NP,base_dir, run_dir, insertion,t1,t2,dur,dur_var,ttl_index);

%% If you want to load LFP:
LFP_Esmolot22 = NP1.getDataLFP(1:384,round(NP1.recordingDuration_ms/2),5000);

%% Get analog data (ECG)

s = NP1.getAnalogData([1 2],0,round(NP1.recordingDuration_ms)-100);

sECG = squeeze(s(2,:));

sDiode = squeeze(s(1,:));
cd(NP1.recordingDir)

save('ECG_Analog.mat',"sECG");
% 
save('AnalogData_Diode.mat','sDiode');

%fb = dwtfilterbank('Wavelet','sym4','SignalLength',numel(s),'Level',3);

%psi = wavelets(fb);
%% Vstim triggers

 [Ttrigger,chNumberT]=NP1.getTrigger();

 Sstart = [Ttrigger{1}];
 Send = [Ttrigger{2}];
 onset = [Ttrigger{3}];
 offset = [Ttrigger{4}];

 %Send = Send(2:end); %For mouse1\Mice_exp_29_11_23\Insertion2\catgt_Mice_exp_29_11_23_2_g0

 ttl_index = 1;

 stimOn = [];
 stimOff = [];
 stimInter = [];
 j =1;


 for i=ttl_index
     stimUp = onset(onset > Sstart(i) & onset < Send(i));
     stimOn = [stimOn stimUp]; %general

     stimDown = offset(offset > Sstart(i) & offset < Send(i));
     stimOff = [stimOff stimDown]; %general

     stimInter = [stimInter mean(stimOff-stimOn)];

     %ttlNum(j) = length(stim(1:2:end)); %sanity check to see how many stimulus presentations there are per round

     j = j+1;
 end

% 
% save('GupTrigger.mat','stimOn')
% save('GdownTrigger.mat','stimOff')
% 
% save('FFFupTrigger.mat','stimOn')
% save('FFFdownTrigger.mat','stimOff')

% GupTrigger = stimOn;
% 
% GdownTrigger = stimOff;
%% triggers Diode
%[stimOn stimOff onSync offSync] = NPdiodeExtract(NP1,1,0,"FFF",1,2); %%%Extract diode using function from diodeextract (special case, no sync)

%%%%Specific case for no sync experiment
fDat=medfilt1(sDiode,600);
% d = designfilt('lowpassiir', 'FilterOrder', 4, ...
%     'HalfPowerFrequency', 120, 'SampleRate', NP1.samplingFrequencyNI);
%%% Apply the filter using filtfilt to preserve phase
% 
% fDat = filtfilt(d, fDat);
stimOn = [];
stimOff = [];
startsFFF = load('Mouse_rightV1_NP_8_7_24_1_g0_tcat.nidq.xd_2_1_0.txt');
endsFFF = load('Mouse_rightV1_NP_8_7_24_1_g0_tcat.nidq.xid_2_1_0.txt');

for i = 1:length(endsFFF)

    signal = fDat(round(startsFFF(i)*NP1.samplingFrequencyNI):round(endsFFF(i)*NP1.samplingFrequencyNI));
    nofiltS = sDiode(round(startsFFF(i)*NP1.samplingFrequencyNI):round(endsFFF(i)*NP1.samplingFrequencyNI));
    stdS = std(signal);

    inSignal = signal*-1;
    inSignal(inSignal>mean(inSignal))= mean(inSignal);
    
    signal(signal>mean(signal)) =mean(signal);

    [pkvalsUp, pklocUp] = findpeaks(signal,'MinPeakProminence',0.8*stdS,'MinPeakDistance',NP1.samplingFrequencyNI*0.05);
    [pkvalsDown, pklocDown] = findpeaks(inSignal,'MinPeakProminence',0.8*stdS,'MinPeakDistance',NP1.samplingFrequencyNI*0.05);
%         figure;plot(signal(2000000:2050000));%hold on;plot(nofiltS(2000000:2050000));
%         xline(pklocUp(pklocUp>2000000 & pklocUp<2050000)-2000000,'k');
%         xline(pklocDown(pklocDown>2000000 & pklocDown<2050000)-2000000,'r');
    stimOn = [stimOn pklocUp];
    stimOff = [stimOff pklocDown];

end

save('FFFupTrigger.mat','stimOn')
save('FFFdownTrigger.mat','stimOff')
%% Plot analog data 

figure
plot(s)
 hold on
% plot(psi(3,:),'r')
plot(-2*circshift(psi(3,:),[0 -38]),'r')
axis tight
legend('QRS Complex','Sym4 Wavelet')
title('Comparison of Sym4 Wavelet and QRS Complex')
hold off
% %% Create filter to extract peaks
%  wt = modwt(s,10);
%  wtrec = zeros(size(wt));
%  inScale=10;
%  wtrec(inScale,:) = wt(inScale,:);
% y = imodwt(wtrec,'sym4');
% sECG
 

%% Extract peaks
% y = abs(y).^2;

Fs = NP1.samplingFrequencyNI; % Example sampling frequency in Hz, adjust this to your signal's sampling frequency
Fn = Fs/2; % Nyquist frequency
Fhp = 0.5; % High-pass filter cutoff frequency in Hz, adjust based on your needs
[b, a] = butter(2, Fhp/Fn, 'high'); % 2nd order Butterworth filter
ecgFiltered = filtfilt(b, a, sECG); % Zero-phase filtering to avoid phase shift

%ecgFiltered=medfilt1(ecgFiltered,500);

[qrspeaks,locs] = findpeaks(ecgFiltered,'MinPeakHeight',mean(ecgFiltered)+1.8*std(ecgFiltered),...
    'MinPeakDistance',300);
% 
% figure()
% plot(ecgFiltered(round(2520*NP1.samplingFrequencyNI):end));yline(mean(ecgFiltered)+std(ecgFiltered));%hold on;plot(normZeroOne(y(1:20000))) 
% hold on
% plot(locs(locs>round(2520*NP1.samplingFrequencyNI))-round(2520*NP1.samplingFrequencyNI),qrspeaks(locs>round(2520*NP1.samplingFrequencyNI)),'ro')
% 
% figure()
% plot(sECG);%hold on;plot(normZeroOne(y(1:20000))) 
% hold on
% plot(locs,qrspeaks,'ro')

locsms = locs/(NP1.samplingFrequencyNI/1000);

%%

%%%% SYNC %%%%%%%%%%%%%%%%%%%%%%


% Search in for neural SQW and NI SQW and load them
cd(NP1.recordingDir)
Neur = readmatrix(dir(fullfile(NP1.recordingDir, '*imec0.ap.xd_384_6_500.txt*')).name);
originSQW = readmatrix(dir(fullfile(NP1.recordingDir, '*nidq.xd_11_0_500.txt*')).name);

locSync = interp1([originSQW(1)-1;originSQW;originSQW(end)+1]'*1000, [Neur(1)-1;Neur;Neur(end)+1]'*1000, locs/(NP1.samplingFrequencyNI/1000), 'linear');



% fileID = fopen("Peak_sync"+".mat",'w');
% fprintf(fileID, '%d\n', locSync);
% fclose(fileID);

save('Peak_nosync.mat','locsms')






%%
figure
plot(y)
hold on
plot(locs,qrspeaks,'ro')
%% Find HR
%Rate

seconds = locs/NP1.samplingFrequencyNI;

r = 1:60:round(max(seconds));
rate = 1:length(r);


for i = 1:length(rate)-1

  rate(i) = sum(seconds>r(i) & seconds<r(i+1));
  
end


%% Input the time of injection

 msdur = NP1.recordingDuration_ms;

 msAfterInject = 741.4395754173818*1000; % got the seconds from the meta of the 4rth recording

 msInject = msdur- msAfterInject; 


%% Build the TIC matrix
p = NP1.convertPhySorting2tIc(NP1.recordingDir); %%Load and save neural data

cd(NP1.recordingDir)
%p = NP1.convertBCSorting2tIc(NP1.recordingDir);


%Select good units
label = string(p.label');

%labelG = label(label == 'good'|label == 'non-somatic');

goodU = p.ic(:,label == 'good');

triggersFFF = load('FFFupTrigger.mat');
triggersFFF = triggersFFF.stimOn;

triggersFFFon = load('Mouse_rightV1_NP_8_7_24_1_g0_tcat.nidq.xd_2_2_0.txt')*1000;

triggersFFFoff = load('Mouse_rightV1_NP_8_7_24_1_g0_tcat.nidq.xid_2_2_0.txt')*1000;

figure;
plot(diff(triggersFFFon))
ylim([0,5000])
bin = 1000;

stimsRange = [4000:4500];

stimDur = mean(triggersFFFoff(stimsRange) - triggersFFFon(stimsRange));

stimInter = mean(triggersFFFon(stimsRange(1:end-1)+1) - triggersFFFoff(stimsRange(1:end-1)));

[M]=BuildBurstMatrix(goodU,round(p.t/bin),round((triggersFFFon(stimsRange)-stimInter/2)/bin)',round((stimDur+stimInter)/bin));


[M2]=BuildBurstMatrix(goodU,round(p.t/bin),round(1/bin)',round(600000/bin));


for u = 1:length(goodU)
    fig = figure;
    imagesc(squeeze(M(:,u,:)));colormap(flipud(gray(64)));
    xline(stimInter/2/bin,'k', LineWidth=1.5)
    xline(stimDur/bin+stimInter/2/bin,'k',LineWidth=1.5)
    ylabel('Trials');xlabel('Time (ms)');
    xticks([0 stimInter/2/bin stimInter/2/bin+stimDur/bin stimInter/bin+stimDur/bin])
    xticklabels(round(xticks-stimInter/2/bin)*bin)
    prettify_plot
    if ~exist(NP1.recordingDir+"\Figs",'dir')
        mkdir Figs
    end
    cd(NP1.recordingDir + "\Figs")
    fig.Position = [1224         259         659         660];
    exportgraphics(fig, sprintf('%s-FFF-ms%d-U%d.pdf',NP1.recordingName,stimDur,u), 'ContentType', 'vector');
    close all

end


fig = figure;
imagesc(squeeze(M2));colormap(flipud(gray(64)));
ylabel('Neurons');xlabel('Time (min)');
xticklabels(round(xticks*bin/1000/60))
prettify_plot
if ~exist(NP1.recordingDir+"\Figs",'dir')
    mkdir Figs
end
cd(NP1.recordingDir + "\Figs")
fig.Position = [1224         259         659         660];
exportgraphics(fig, sprintf('%s-FFF-ms%d-allUnits-10min.pdf',NP1.recordingName,stimDur), 'ContentType', 'vector');
close all


% 
% spike_times = double(readNPY('spike_times_original.npy'))*(NP1.samplingFrequency)/(30000);
% 
% writeNPY(spike_times,'spike_times.npy');

%*30000/NP1.samplingFrequencyAP


save('SpikesPerNeuron.mat','M','-v7.3') %%Multiply times as Mark says.

%%
bin=2;
win= 160;
start = 40;
%% depth
cluster_info = readtable(string(NP1.recordingDir) + "\cluster_info.tsv",  "FileType","text",'Delimiter', '\t');
GoodU_or = cluster_info.cluster_id(cluster_info.group=="good");
GoodU_orDepth = cluster_info.depth(cluster_info.group=="good");

verticalDepth = sin(deg2rad(72.5))*(2832 - GoodU_orDepth);


cd(NP1.recordingDir)
%save('unitDepth',"verticalDepth");

%% times before esmolol injection
hb = (locs/(NP1.samplingFrequencyNI/1000))';

hb1 = hb(hb < msInject);

%% Build the burst matrix Before esmolol 



[nTrials,nNeurons,nTimes]=size(M);
%%
figure;
subplot(2,1,1)
%Plot average accros trials where Units are the y axis. 
%imagesc((1:nTimes)*bin,1:nNeurons,normalize(squeeze(mean(M,1)), 2, 'range'));colormap(flipud(gray(256)));ylabel('Neurons');
imagesc((1:nTimes)*bin,1:nNeurons,squeeze(mean(M,1)));colormap(flipud(gray(256)));ylabel('Neurons');
%title(title_stim); 
xlim([1 nTimes*bin]);
cb=colorbar;set(cb,'Position',[0.91 0.75 0.013 0.17]);ylabel(cb,'Spk. Rate');
line([100 100],ylim,'color','r'); 
%yticklabels(string(goodU(1,:))); 
%yticks(1:nNeurons); 
set(gca,'YTick',[]);

%% times after esmolol injection
hb2 = hb(hb > msInject);

[M]=BuildBurstMatrix(goodU,round(p.t/bin),round(hb2/bin)-100,round(win/bin)); %leave noise units behind

[nTrials,nNeurons,nTimes]=size(M);
%%
figure;
subplot(2,1,1)
%Plot average accros trials where Units are the y axis. 
imagesc((1:nTimes)*bin,1:nNeurons,normalize(squeeze(mean(M,1)), 2, 'norm', Inf));colormap(flipud(gray(256)));ylabel('Neurons');
%imagesc((1:nTimes)*bin,1:nNeurons,squeeze(mean(M,1)));colormap(flipud(gray(256)));ylabel('Neurons');%
%without normalization
%title(title_stim); 
xlim([1 nTimes*bin]);
cb=colorbar;set(cb,'Position',[0.91 0.75 0.013 0.17]);ylabel(cb,'Spk. Rate');
line([100 100],ylim,'color','r'); 
%yticklabels(string(goodU(1,:))); 
%yticks(1:nNeurons); 
set(gca,'YTick',[]);
%% Run for M before and after esmolol inject (hb1 and hb2)
hbRaw = NP1.getAnalogData(1,hb2-start,win);

[~,hbs,samples] = size(hbRaw);

hbRawM = squeeze(mean(hbRaw, 2))';


x = 1:length(hbRawM);
dx = linspace(min(x),max(x),round(samples/(NP1.samplingFrequencyNI/1000)));

dy = interp1(x, hbRawM, dx); 

msecondsX = (dx/(NP1.samplingFrequencyNI/1000));

%% 
subplot(2,1,2)
plot(msecondsX, dy, '-r');
ylabel('Volts');xlabel('Time [ms]');
xlim([0 max(msecondsX)]);
xline(100,'color','b'); 
xticklabels(0:20:length(dx));

%%
hold on
NormSpikes = normalize(squeeze(mean(M,1)), 2, 'range');
sumNS = sum(NormSpikes,1);
sumred = movsum(sumNS,5);
sumredF = sumred(3:5:end)/max(sumred(3:5:end));
hold on
%sumredF = rand(length(sumredF),1);
b1 = bar(1:10:round(max(msecondsX)),sumredF);
b1.FaceAlpha = 0.5;

%% plot ind neurons.

pathfig = '\\132.66.45.127\data\Large_scale_mapping_NP\Mice_experiments\Esmolo_test_2\Insertion2\\Figures\after_esmolol';

cd(pathfig)

rands = 1000; %Randomization permutations

ti = 'After esmolol';

%in = [3, 92, 5, 6,10,13,15,18,21,26,28,29,30,34,41,45,49,57,60,63,65,67,68,74,76,79,83,86,88,89,90,91,92,94,95,96,99,]; 
for i = 1:length(goodU) %1:105

    Mu = M(:,i,:); 
    [nTrials,nNeurons,nTimes]=size(Mu);
    title_stim = sprintf('Unit-%d-channel-#%s-%s', i, string(goodU(1,i)),ti);


    shuf = zeros(rands,size(Mu,3));
    for s = 1:rands

        M_shufS=Mu(:,:,randperm(size(Mu,3)));
        shuf(s,:) = sum(squeeze(M_shufS));

    end
    
    sumHB = movsum(squeeze(Mu), [9 0]);
    
    sumS = sumHB(10:10:end,:);
    %figure('visible','off'); 
    figure;
    subplot(3,1,1:2);
    imagesc((1:nTimes)*bin,1:nTrials,sumS);colormap(flipud(gray(256)));ylabel('# heartbeats');xlabel('Time [ms]');
    title(title_stim); 
    cb=colorbar;set(cb,'Position',[0.91 0.75 0.013 0.17]);ylabel(cb,'Spk. Rate');
    line([start start],ylim,'color','r'); 
    x0=10;
    y0=10;
    width=1000;
    height=hbs*2;
    set(gcf,'position',[x0,y0,width,height])
    
    subplot(3,1,3)
    plot(msecondsX, dy, '-r', LineWidth=1.5);
    ylabel('Volts');xlabel('Time [ms]');
    xlim([0 max(msecondsX)]);
    xline(start,'color','b');
    %xticklabels(0:20:length(dx));

    hold on

    %Shuffled data:

    %For testing
    B = reshape(shuf, rands, 5, []); % Reshape it to a 100x5x16 matrix
    C = sum(B, 2); % Sum along the 2nd dimension
    S = squeeze(C); % Remove the singleton dimensions


    %For plotting
    SmSpikes =mean(shuf,1);
    [shufC Sedges] = histcounts(SmSpikes);
    SsumSpikes = movsum(SmSpikes, [4 0]);
    Ssumred = SsumSpikes(5:5:end);

    %Real data:
    mSpikes = sum(squeeze(Mu));
    [C edges] = histcounts(mSpikes);
    sumSpikes = movsum(mSpikes, [4 0]);
    sumred = sumSpikes(5:5:end);

    %Statistical test:
    pvals = zeros(1,length(sumred));
    h = zeros(1,length(sumred));
    for p = 1:length(sumred)
        %[h(p) pvals(p)] = ttest(S(:,p)',sumred(p),'Tail','right');
        
        pvalsH(p) = (sum(S(:,p)' >= sumred(p)) +1)/(rands+1);
        pvalsL(p) = (sum(S(:,p)' <= sumred(p)) +1)/(rands+1);

        if (pvalsH(p) < 0.05) || (pvalsL(p) < 0.05) 
            h(p) = 1;
        end

%         figure;
%         histogram(S(:,p)');
%         xline(sumred(p))
    end


    x = 10:10:round(max(msecondsX));
    yd = sumred/max(sumred)-0.05;

    b1 = bar(x,yd,'FaceColor',[0 .5 .5],'EdgeColor','white');
    b1.FaceAlpha = 0.5;
    

    for p = 1:length(pvals)

        if h(p) ==1
            
            text(x(p)-(x(2)-x(1))/8,yd(p), '*','FontSize',20)

        end
    end
    
    ys = Ssumred/max(sumred)-0.05;

    b2 = bar(x,ys,'FaceColor','white','EdgeColor','black','LineStyle','--');
    b2.FaceAlpha = 0;
    b2.EdgeAlpha = 0.5;




    %export_fig test.png -transparent
    print(title_stim,'-depsc2');
    print(title_stim,'-dpng');
    %xticklabels(0:20:length(dx));
    % yticks(1:nNeurons); 
    close all

end
%%
winMin = 200;
NeuronID = 30;
j = 1;
for i = 1:10000:round(NP1.recordingDuration_ms)-100
        
    hbmin = hb((i<hb) & (hb<i+60000));

    
    [Mm]=BuildBurstMatrix(goodU(:,NeuronID),round(p.t/bin),round(hbmin/bin)-100,round(winMin/bin));

    mSpikes = sum(squeeze(Mm));

    sumSpikes = movsum(mSpikes, [2 0]);

    sumred(:,j) = sumSpikes(3:3:end);
    
    j = j+1;


end
%%
[x, y] = size(sumred);

imagesc(sumred);
yline(16.5, 'r', 'Heart beat', LineWidth=1.5)
xline(msInject/10000, 'r','Esmolol injection', LineWidth=1.5);
yticklabels(30:30:200);
xticks(6:30:y)
xticklabels(1:5:22)
xlabel('Time (min)')
ylabel('Time (ms)')


%%
A = readNPY(string(NP1.recordingDir) + "\amplitudes.npy");


cluster_info = readtable(string(NP1.recordingDir) + "\cluster_info.csv"); %remember to convert tsv to csv file

spikeNames = readNPY(string(NP1.recordingDir) + "\spike_clusters.npy");

spikeTimes = readNPY(string(NP1.recordingDir) + "\spike_times.npy")/(NP1.samplingFrequency/1000);

g = cluster_info(matches(cluster_info.group,'good'),:);

good_ids = g{:,1};

good_ids(92)
%%

Amps = A(spikeNames == 219);
Atimes = cast(spikeTimes(spikeNames == 219),'double'); %transform to mins

plot(Atimes,Amps);
xline(msInject, 'r','Esmolol injection', LineWidth=1.5);

%%
hold on
%Get coefficients of a line fit through the data.
coefficients = polyfit(Atimes, Amps, 1);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(Atimes), max(Atimes), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);
% Plot everything.
plot(Atimes, Amps, 'b.', 'MarkerSize', 15); % Plot training data.
hold on; % Set hold on so the next plot does not blow away the one we just drew.
plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot fitted line.
grid on;

%%
%Plot raw analog data against 

tNDQ = load('\\132.66.45.127\data\Large_scale_mapping_NP\Mice_experiments\Esmolo_test_2\Insertion2\catgt_Esmolo_test_2_2_g0\Esmolo_test_2_2_g0_tcat.nidq.xd_16_0_500.txt');

tNP = load('\\132.66.45.127\data\Large_scale_mapping_NP\Mice_experiments\Esmolo_test_2\Insertion2\catgt_Esmolo_test_2_2_g0\Esmolo_test_2_2_g0_tcat.imec0.ap.xd_384_6_500.txt');

diff = abs(tNP-tNDQ);


plot(1:length(tNP),diff*1000)
ylabel('miliseconds');
xlabel('square-wave event');
%%
timeSeriesViewer(NP1)


