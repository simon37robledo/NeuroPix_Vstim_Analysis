%Sleep analysis
basic_path = '\\132.66.45.127\data\Large_scale_mapping_NP\\Awake_exp\SA6';
exp = 'SA6_experiment_4_9_23';
in =4;

path = convertStringsToChars(basic_path+"\"+exp+"\Insertion"+in+"\catgt_"+exp+"_"+in+"_g0");

NP = NPAPRecording(path);
%%
timeSeriesViewer(NP)
%%

p = NP.convertPhySorting2tIc(NP.recordingDir);

%Select good units
label = string(p.label');

goodU = p.ic(:,label == 'good');

bin =20;
%M = squeeze(BuildBurstMatrix(goodU,round(p.t/bin),round(6802000/bin),round(100000/bin)));

M = squeeze(BuildBurstMatrix(goodU,round(p.t/bin),round(6846000/bin),round(5000/bin)));

cluster_info = readtable(string(NP.recordingDir) + "\cluster_info.tsv",  "FileType","text",'Delimiter', '\t');

cluster_group = readtable(string(NP.recordingDir) + "\cluster_group.tsv",  "FileType","text",'Delimiter', '\t');

GoodU_or = cluster_info.cluster_id(cluster_info.group=="good");

goodU_chans = cluster_info.ch(cluster_info.group=="good");

idx = find(goodU_chans>240 & goodU_chans < 260);

[nN, nB] = size(M);

ch = cluster_info.ch(cluster_info.cluster_id == 304);

imagesc()

%%

LFP = NP.getData(1,5,300);

figure;
plot(squeeze(LFP))

LFP2 = NP.getDataLFP(1,5,300);

figure;
plot(squeeze(LFP2))
%%
figure(7);
limits = axis;
hold on
rectangle('Position', [44000,limits(3) , 5000, limits(4)*3], 'FaceColor',[0 0 0 0.2]);

%%
figure
%imagesc(M)

imagesc((1:nB)*bin,1:nN,M);colormap(flipud(gray(64)));ylabel('Neurons');xlabel('Time (ms)');caxis([0 1])

%%
d = NP.getData(250, 4500000,4500000+30000);

d =  NP.getData(250:255, 4500000,4500000+30000);



d = squeeze(d);
%plot(d')

l = NP.getDataLFP(250,0,30000);

figure()

plot(squeeze(d(1,:,:)));
xticklabels([0:5:30])

%plot(squeeze(d(1,1,1:10000)));
% 
%timeSeriesViewer(NP)

%%
d =  NP.getData(250:255, 0,30000);


%%

SA=sleepAnalysis('\\132.66.45.127\data\Large_scale_mapping_NP\\Awake_exp\SA6\SleepSA.xlsx');
SA.setCurrentRecording('Animal=SA6,recNames=Insertion4');
DB=SA.getDelta2BetaRatio;
SA.plotDelta2BetaRatio;

%%

startTime = round(NP.recordingDuration_ms - 35*60*1000);

window = 34*60*1000;

DB2 = SA.getDelta2BetaAC('tStart',1000*60*60,'win', NP.recordingDuration_ms-1000*60*60);

SA.plotDelta2BetaAC;

SA.plotDelta2BetaSlidingAC



%DB2 = SA.getDelta2BetaAC;



