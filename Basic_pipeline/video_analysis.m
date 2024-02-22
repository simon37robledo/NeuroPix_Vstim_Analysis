%% Process video files

%initialize class
path = '\\132.66.45.127\data\Large_scale_mapping_NP\\Awake_exp\SA5\SA5_Experiment_10_05_23\Insertion1\catgt_SA5_Experiment_10_05_23_1_g0';


NP = NPAPRecording(path);


[Ttrigger,chNumberT]=NP.getTrigger();

%%
pathV = '\\132.66.45.127\data\Large_scale_mapping_NP\\Awake_exp\SA5\SA5_Experiment_10_05_23\Insertion1\Instertion1_Camera1_20230510-163510.csv';
VideoTS = readtable(pathV); 


TTL_TS = [Ttrigger{4}];

TTL_TS = TTL_TS(1:2:end); %Get only "on"s

idx_StartEnd = find(TTL_TS(2:end) - TTL_TS(1:end-1) > 100); %Get index of 1 second gaps that indicate start and end of recording

plot(diff(TTL_TS)) %Plot to check if there are indeed two gaps

TTL_TSf = TTL_TS(idx_StartEnd(1)+1:idx_StartEnd(2)); % Select TTLs between these two gaps

%Check if number of TTLs concide with frams in VideoTS


check = diff(TTL_TSf) - diff(VideoTS{:,:})'.*1000; %Difference of difference between frames. It should be close to 0.0 ms



