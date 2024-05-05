
%1. %% PREPROCESSING %%%%%%%%%%%%%
%a) Create appropiate folders in local VSComp (PV27, experiment, insertion)
%d) Transfer neural and vstim data to KS PC. 
%b) Create appropiate variables (path, digital channels) from excel file 
%c) Run Preprocessing script, spikeGLX cat_GT with filters

%2. %% SPIKE SORTING %%%%%%%%%%%%%
%a) Run KS3-4 from script
%b) Upload results to SIL3

%3. %% POSTPROCESSING %%%%%%%%%%%%
%a) Run spike duplicates removal

%% Run on old experiments

cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

for i =1:size(data,1)
    path = convertStringsToChars(string(data.Base_path(i))+filesep+string(data.Exp_name(i))+filesep+"Insertion"+string(data.Insertion(i))...
        +filesep+"catgt_"+string(data.Exp_name(i))+"_"+string(data.Insertion(i))+"_g0");
    cd(path)
    NP = NPAPRecording(path);
    
    if ~exist(path+"\oldRez",'dir')
        mkdir oldRez
        movefile rez.mat oldRez
        cd(path+"\oldRez")
    else
        rez = load('rez.mat');rez = rez.rez;
    end

    rez = remove_ks2_duplicate_spikes(rez);
    rez = find_merges(rez,1);
    save('rez.mat','rez')

    rezToPhy2(rez,NP.recordingDir)


end

%%
%b) Make pc-features
%c) Compute spike positions
%d) Run BC?


