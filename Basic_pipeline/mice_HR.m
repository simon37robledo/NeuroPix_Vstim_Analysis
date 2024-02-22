%Mice HR analysis

%Load recording and create neuropixel class

base_dir = "\\132.66.45.127\data\Large_scale_mapping_NP\Mice_experiments\Esmolo_test_2"; %change to NP.recordingDir
run_dir = "\\132.66.45.127\data\Large_scale_mapping_NP\SpikeGLX\"; %Folder that has catgt and tprime subfolders
insertion = "2";
t1 = 0.8;
t2 = 0.9;
dur = 12;
dur_var = 5;
ttl_index =-1;

recPath1 = '\\132.66.45.127\data\Large_scale_mapping_NP\Mice_experiments\Esmolo_test_2\Insertion2\Filter_catgt_Esmolo_test_2_2_g0'; 
NP1= NPAPRecording(recPath1);


%[astim, astimon, astimoff, astimdur,astiminter] = time_diode_synced(NP,base_dir, run_dir, insertion,t1,t2,dur,dur_var,ttl_index);

%% If you want to load LFP:
LFP_Esmolot22 = NP1.getDataLFP(1:384,round(NP1.recordingDuration_ms/2),5000);

%% Get analog data (ECG)

s = NP1.getAnalogData(1,0,round(NP1.recordingDuration_ms)-100);

s = squeeze(s);

%fb = dwtfilterbank('Wavelet','sym4','SignalLength',numel(s),'Level',3);

%psi = wavelets(fb);

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
%% Create filter to extract peaks
 wt = modwt(s,10);
 wtrec = zeros(size(wt));
 inScale=10;
 wtrec(inScale,:) = wt(inScale,:);
y = imodwt(wtrec,'sym4');
plot(normZeroOne(s));hold on;plot(normZeroOne(y))     

%% Extract peaks
% y = abs(y).^2;
[qrspeaks,locs] = findpeaks(s,'MinPeakHeight',0.76,...
    'MinPeakDistance',500);
%%
figure
plot(s)
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
p = NP1.convertPhySorting2tIc(NP1.recordingDir);

%Select good units
label = string(p.label');

goodU = p.ic(:,label == 'good');

bin=2;
win= 160;
start = 40;

%% times before esmolol injection
hb = (locs/(NP1.samplingFrequencyNI/1000))';

hb1 = hb(hb < msInject);

%% Build the burst matrix Before esmolol 

[M]=BuildBurstMatrix(goodU,round(p.t/bin),round(hb1/bin)-start,round(win/bin));

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


