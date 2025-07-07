%%%%Preprocess function
cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

%%

%%% Call preprocessing:

for ex =[79:83]
    experiment = data(ex,:);
    FRunCatGTnTprime(experiment);

end

function FRunCatGTnTprime(data)


run_dir1 = "\\sil1\data\Large_scale_mapping_NP\SpikeGLX\V4_catGT";
run_dir2 = "Y:\Large_scale_mapping_NP\SpikeGLX\V4_catGT";

% run_dir1 = "Y:\data\Large_scale_mapping_NP\SpikeGLX\V4_catGT";
% run_dir1 = "\\169.254.69.105\data\Large_scale_mapping_NP\SpikeGLX\V4_catGT";
% run_dir1 = "\\169.254.248.184\data\Large_scale_mapping_NP\SpikeGLX\V4_catGT";

%% %% Phy commmand: phy template-gui params.py

% basic_pathPV102 = '\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\PV102';
% expPV102 = 'PV102_experiment_18_7_23';
% 
% basic_pathPV103 = '\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\PV103';
% expPV103 = 'PV103_Experiment_12_6_23';
% 
% basic_pathPV67 = '\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\PV67';
% expPV67= 'PV67_experiment_5_7_23';
% 
% basic_pathPV27 = '\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\PV27';
% expPV27 = 'PV27_Experiment_25_6_23';
% 
% basic_pathPV139 = '\\sil3\data\Large_scale_mapping_NP\lizards\PV139';
% expPV139 = 'PV139_Experiment_6_2_24';
% 
% basic_pathPV59 = '\\sil3\data\Large_scale_mapping_NP\lizards\PV59';
% expPV59 = 'PV59_Experiment_20_2_24';
% 
% basic_pathPV32 = '\\sil3\data\Large_scale_mapping_NP\lizards\PV32';
% expPV32 = 'PV32_Experiment_18_3_24';
% 
% basic_pathPV32 = '\\sil3\data\Large_scale_mapping_NP\lizards\PV32';
% expPV32 = 'PV32_Experiment_18_3_24';
% 
% basic_pathPV152 = '\\sil3\data\Large_scale_mapping_NP\lizards\PV152';
% expPV152 = 'PV152_Experiment_11_7_24';
% 
% basic_pathPV43 = '\\sil3\data\Large_scale_mapping_NP\lizards\PV43';
% expPV43 = 'PV43_Experiment_24_7_24';
% 
% 
% basic_pathPV104 = '\\sil3\data\Large_scale_mapping_NP\lizards\PV104';
% expPV104= 'PV104_Experiment_5_8_24';
% 
% basic_pathPV35 = '\\sil3\data\Large_scale_mapping_NP\lizards\PV35';
% expPV35 = 'PV35_Experiment_18_8_24';
% 
% basic_pathSA8 = '\\sil3\data\Large_scale_mapping_NP\lizards\SA8';
% expSA8 = 'SA8_Experiment_7_11_24';
% 
% basic_pathPV97= '\\sil3\data\Large_scale_mapping_NP\lizards\PV97';
% expPV97= 'PV97_Experiment_21_01_25';

basic_path = data.Base_path;
expPath = data.Exp_name;


%% Basic variables
base_dir = string(basic_path)+"\"+string(expPath);
%Folder that has catgt and tprime subfolders
insertion = data.Insertion;
%fileName = "PV97_Experiment";
runs = "mult";
dig_CH = string(data.Digital_channel);
concat = 1;
syncChan = string(data.Sync_bit);


%% Command line excecution

try
    cd(base_dir)
catch
    originP = cell2mat(extractBetween(base_dir,"\\","\Large_scale"));
    if strcmp(originP,'sil3\data')
        base_dir = replaceBetween(base_dir,"","\Large_scale","W:");
    else
        base_dir = replaceBetween(base_dir,"","\Large_scale","Y:");
    end
    cd(base_dir)
end


%0.0 Set directory to catGT directory
try
    cd(run_dir1+"\CatGT-win")
catch
    cd(run_dir2+"\CatGT-win")
end

%0.1 Find experiment file name from base_dir
out=regexp(base_dir,'\','split');
exp = string(out(end));
insertion = string(insertion);

%0.2 Find number of stimulus in insertion file
if insertion ~= "-1"
    file = dir (base_dir + "\Insertion" + insertion);
else
    file = dir (base_dir);
end


filenames = {file.name};
num = string(sum( ~cellfun(@isempty, strfind(filenames, expPath)))-1); %Normal is -1

%%
if concat ==1

%1. Create command line filter + concatenation ap
if insertion ~= "-1"
    

    cmndAP = "CatGT -dir="+ base_dir + "\Insertion" + insertion + " -run=" + exp +"_" + insertion ...
        + " -g=0," + num + " -t=0,0 -prb_fld -t_miss_ok " + " -ap" + " -prb=0 -zerofillmax=0"...
        + " -dest="...
        + base_dir + "\Insertion" + insertion;

    status = 1;

    %1.1. execute command
    status = system(cmndAP);

    while status == 1
        %waiting for command to complete
    end

    disp("AP extraction and concatenation completed")

else
    cmndAP = "CatGT -dir="+ base_dir + " -run=" + exp ...
        + " -g=0," + num + " -t=0,0 -prb_fld -t_miss_ok " + " -ap" + " -prb=0 -zerofillmax=0"...  %+ " -apfilter=butter,12,300,9000 -gfix=0.40,0.10,0.02 -loccar=4,3
        + " -dest="...
        + base_dir;

    status = 1;

    %1.1. execute command
    status = system(cmndAP);

    while status == 1
        %waiting for command to complete
    end

    disp("AP extraction completed")


end
else

    cd(base_dir + "\Insertion" + insertion)

    mkdir("catgt_" + exp + "_" + insertion + "_g0" + "\" + exp + "_" + insertion);

    source = base_dir + "\Insertion" + insertion + "\" + exp + "_" + insertion + "_g0\" + exp + "_" + insertion + "_g0_imec0";

    dest = base_dir + "\Insertion" + insertion + "\" + "catgt_" + exp + "_" + insertion + "_g0" + "\" + exp + "_" + insertion + "_g0_tcat.imec0.";

    cd(source)
    movefile(exp + "_" + insertion + "_g0_t0.imec0.ap.bin",dest+"ap.bin")
    movefile(exp + "_" + insertion + "_g0_t0.imec0.ap.meta",dest + "ap.meta")
    movefile(exp + "_" + insertion + "_g0_t0.imec0.lf.bin",dest + "lf.bin")
    movefile(exp + "_" + insertion + "_g0_t0.imec0.lf.meta",dest + "lf.meta")

    movefile(dest + "ap.bin",source+"\" + exp + "_" + insertion + "_g0_t0.imec0.ap.bin")
    movefile(dest + "ap.meta", source+"\"+exp + "_" + insertion + "_g0_t0.imec0.ap.meta")
    movefile(dest + "lf.bin",source+"\"+exp + "_" + insertion + "_g0_t0.imec0.lf.bin")
    movefile(dest + "lf.meta",source+"\"+exp + "_" + insertion + "_g0_t0.imec0.lf.meta")
    
    disp("binary and meta files moved to catgt folder")

end



%2. Create command line filter + concatenation lf
if concat == 1

if insertion ~= "-1"

    cmndLF = "CatGT -dir="+ base_dir + "\Insertion" + insertion + " -run=" + exp +"_" + insertion ...
        + " -g=0," + num + " -t=0,0 -prb_fld -t_miss_ok " + " -lf" + " -prb=0 -zerofillmax=0"...
        + " -dest="...
        + base_dir + "\Insertion" + insertion;

    status = 1;

    %2.1. execute command
    status = system(cmndLF);

    while status == 1
        %waiting for command to complete
    end

    disp("LF extraction and concatenation completed")

else
    cmndLF = "CatGT -dir="+ base_dir + " -run=" + exp ...
        + " -g=0," + num + " -t=0,0 -prb_fld -t_miss_ok " + " -lf" + " -prb=0 -zerofillmax=0"...
        + " -dest="...
        + base_dir;

    status = 1;

    %2.1. execute command
    status = system(cmndLF);

    while status == 1
        %waiting for command to complete
    end

    disp("LF extraction and concatenation completed")

end
end


%3. Create command line for concatenation nidaq

if insertion == "-1" & runs == "single"

    disp("Single run, no concatenation")

elseif insertion == "-1" & runs ~= "single"
    cmndNIc =  "CatGT -dir="+ base_dir + " -run=" + exp ...
        + " -g=0," + num + " -t=0,0 -prb_fld -t_miss_ok " + " -ni" + " -prb=0 -zerofillmax=0"...
        + " -dest="...
        + base_dir;

    status = 1;

    %3.1. execute command
    status = system(cmndNIc);

    while status == 1
        %waiting for command to complete
    end

    disp("NI concatenation completed")
else
    cmndNIc = "CatGT -dir="+ base_dir + "\Insertion" + insertion + " -run=" + exp +"_" + insertion ...
        + " -g=0," + num + " -t=0,0 -prb_fld -t_miss_ok " + " -ni" + " -prb=0 -zerofillmax=0"...
        + " -dest="...
        + base_dir + "\Insertion" + insertion;

    status = 1;

    %3.1. execute command
    status = system(cmndNIc);

    while status == 1
        %waiting for command to complete
    end

    disp("NI concatenation completed")


end


% If the number of stimulations is 1, then cut, paste and rename the
% nidq file.
if num == "0" && insertion ~="-1"
    source = base_dir + "\Insertion" + insertion + "\" + exp + "_" + insertion + "_g0";
    dest = base_dir + "\Insertion" + insertion + "\" + "catgt_" + exp + "_" + insertion + "_g0" + "\" + exp + "_" + insertion + "_g0_tcat.nidq.bin";
    cd(source)

    copyfile(exp+"_"+ insertion + "_g0_t0.nidq.bin",dest);

    %     elseif insertion == -1
    %         source = base_dir + "\" + exp + "_g0";
    %         dest = base_dir + "\" + "catgt_" + exp + "_g0" + "\" + exp + "_g0_tcat.nidq.bin";
    %         cd(source);
    %         copyfile(exp+"_"+ insertion + "_g0_t0.nidq.bin",dest);


end


try
    cd(run_dir1+"\CatGT-win")
catch
    cd(run_dir2+"\CatGT-win")
end

%4. Create command line for extracting nidaq

if insertion ~= "-1"

    cmndNIe = "CatGT -dir="+ base_dir + "\Insertion" + insertion + " -run=catgt_" + exp +"_" + insertion ...
        + " -g=0 -ni -prb=0 -t=cat -no_tshift "...
        + "-xd=0,0,-1,0,0 "...
        + "-xd=0,0,-1,1,0 "...
        + "-xd=0,0,-1,2,0 "...
        + "-xd=0,0,-1,3,0,1 "...
        + "-xd=0,0,-1,5,0 "...
        + "-xid=0,0,-1,1,0 "...
        + "-xid=0,0,-1,2,0 "...
        + "-xid=0,0,-1,3,0,1 "... %% The signal is req1uired to stay low for at least 1 sample (When default, 5, missing down triggers
        + "-xid=0,0,-1,5,0 "... %%cAMERA TTLs
        + "-dest="...
        + base_dir + "\Insertion" + insertion;

    status = 1;

    %4.1. execute command
    status = system(cmndNIe);

    while status == 1
        %waiting for command to complete
    end

    disp("NI extraction completed")
else
    cmndNIe = "CatGT -dir="+ base_dir + " -run=catgt_" + exp ...
        + " -g=0 -ni -prb=0 -t=cat -no_tshift "...
        + "-xd=0,0,-1,0,0 "...
        + "-xd=0,0,-1,1,0 "...
        + "-xd=0,0,-1,2,0 "...
        + "-xd=0,0,-1,3,0,1 "...
        + "-xd=0,0,-1,5,0 "...
        + "-xid=0,0,-1,1,0 "...
        + "-xid=0,0,-1,2,0 "...
        + "-xid=0,0,-1,3,0,1 "...
        + "-xid=0,0,-1,5,0 "...
        + "-dest="...
        + base_dir;

    status = 1;

    %4.1. execute command
    status = system(cmndNIe);

    while status == 1
        %waiting for command to complete
    end

    disp("NI extraction completed")

end


%5.0. Create sync folder

if insertion ~= "-1"

    path = base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0";
else
    path =  base_dir + "\catgt_" + exp + "_g0";
end

cd(path)

if ~exist(path+"\sync_events",'dir')
    mkdir sync_events
end

%Go to TPrime dir


try
    cd(run_dir1 + "\TPrime-win")
catch
    cd(run_dir2 + "\TPrime-win")
end




%5.1. Create TPrime command

if insertion ~= "-1"
    cmndTPrime = "TPrime -syncperiod=1.0 -tostream=" + base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + exp + "_" + insertion + "_g0_tcat.imec0.ap.xd_384_"+"6"+"_500.txt "...
        + "-fromstream=0,"+ base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + exp + "_" + insertion + "_g0_tcat.nidq.xd_"+ dig_CH+"_"+syncChan+"_500.txt "...
        + "-events=0,"+ base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + exp + "_" + insertion + "_g0_tcat.nidq.xd_"+ dig_CH+"_1_0.txt,"...
        + base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + "sync_events\out_1.txt "... %stimulus onset
        + "-events=0,"+ base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + exp + "_" + insertion + "_g0_tcat.nidq.xd_"+ dig_CH+"_2_0.txt,"...
        + base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + "sync_events\out_2.txt "... %stimulus trial onset
        + "-events=0,"+ base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + exp + "_" + insertion + "_g0_tcat.nidq.xd_"+ dig_CH+"_3_0.txt,"...
        + base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + "sync_events\out_3.txt "... %digital frame onset
        + "-events=0,"+ base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + exp + "_" + insertion + "_g0_tcat.nidq.xd_"+ dig_CH+"_5_0.txt,"...
        + base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + "sync_events\out_5.txt "... %camera frame onset
        + "-events=0,"+ base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + exp + "_" + insertion + "_g0_tcat.nidq.xid_"+ dig_CH+"_1_0.txt,"...
        + base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + "sync_events\out_1inv.txt "... %stimulus offset
        + "-events=0,"+ base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + exp + "_" + insertion + "_g0_tcat.nidq.xid_"+ dig_CH+"_2_0.txt,"...
        + base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + "sync_events\out_2inv.txt "... %stimulus trial offset
        + "-events=0,"+ base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + exp + "_" + insertion + "_g0_tcat.nidq.xid_"+ dig_CH+"_3_0.txt,"...
        + base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + "sync_events\out_3inv.txt "... %digital frame offset
        + "-events=0,"+ base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + exp + "_" + insertion + "_g0_tcat.nidq.xid_"+ dig_CH+"_5_0.txt,"...
        + base_dir + "\Insertion" + insertion + "\catgt_" + exp + "_" + insertion + "_g0\" + "sync_events\out_5inv.txt"; %camera frame offset
else
    cmndTPrime = "TPrime -syncperiod=1.0 -tostream=" + base_dir + "\catgt_" + exp + "_g0\" + exp + "_g0_tcat.imec0.ap.xd_384_6_500.txt "...
        + "-fromstream=0,"+ base_dir  + "\catgt_" + exp + "_g0\" + exp + "_g0_tcat.nidq.xd_"+ dig_CH+"_0_500.txt "...
        + "-events=0,"+ base_dir + "\catgt_" + exp + "_g0\" + exp + "_g0_tcat.nidq.xd_"+ dig_CH+"_1_0.txt,"...
        + base_dir +  "\catgt_" + exp + "_g0\" + "sync_events\out_1.txt "...
        + "-events=0,"+ base_dir + "\catgt_" + exp + "_g0\" + exp + "_g0_tcat.nidq.xd_"+ dig_CH+"_2_0.txt,"...
        + base_dir + "\catgt_" + exp + "_g0\" + "sync_events\out_2.txt "...
        + "-events=0,"+ base_dir + "\catgt_" + exp + "_g0\" + exp + "_g0_tcat.nidq.xd_"+ dig_CH+"_3_0.txt,"...
        + base_dir +  "\catgt_" + exp + "_g0\" + "sync_events\out_3.txt "...
        + "-events=0,"+ base_dir + "\catgt_" + exp + "_g0\" + exp + "_g0_tcat.nidq.xd_"+ dig_CH+"_5_0.txt,"...
        + base_dir +  "\catgt_" + exp + "_g0\" + "sync_events\out_5.txt "...
        + "-events=0,"+ base_dir + "\catgt_" + exp +  "_g0\" + exp  + "_g0_tcat.nidq.xid_"+ dig_CH+"_1_0.txt,"...
        + base_dir  + "\catgt_" + exp + "_g0\" + "sync_events\out_1inv.txt "...
        + "-events=0,"+ base_dir + "\catgt_" + exp  + "_g0\" + exp + "_g0_tcat.nidq.xid_"+ dig_CH+"_2_0.txt,"...
        + base_dir + "\catgt_" + exp +  "_g0\" + "sync_events\out_2inv.txt "...
        + "-events=0,"+ base_dir + "\catgt_" + exp  + "_g0\" + exp + "_g0_tcat.nidq.xid_"+ dig_CH+"_3_0.txt,"...
        + base_dir + "\catgt_" + exp +  "_g0\" + "sync_events\out_3inv.txt "...
        + "-events=0,"+ base_dir + "\catgt_" + exp  + "_g0\" + exp + "_g0_tcat.nidq.xid_"+ dig_CH+"_5_0.txt,"...
        + base_dir + "\catgt_" + exp +  "_g0\" + "sync_events\out_5inv.txt";
end




%5.2. execute command
status = system(cmndTPrime);

while status == 1
    %waiting for command to complete
end

disp("Event syncing completed, look for files in sync_events folder")

end