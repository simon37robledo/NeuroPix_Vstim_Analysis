%%Run all stimuli scripts
SDGrecordingsA = [8:14,40:43,49:54]; %anesthetized 

GoodRecordingsPV =[8:21,40:43,49:54];

%Visual stim computer
%files = dir('C:\Users\MarkS9\Documents\GitHub\NeuroPix_Vstim_Analysis\Basic_pipeline\StimulusAnalysis\*.m'); 

%windows server
files = dir('D:\Mark_S13\Documents\GitHub\NeuroPix_Vstim_Analysis\Basic_pipeline\StimulusAnalysis\*.m');

filenames = {files.name};
fileIndex = find(~contains(filenames,'runAllstims') & (contains(filenames,'Moving_ball')|contains(filenames,'Rectangle_grid')|contains(filenames,'stimGratings')));% Avoid infinite loop and select stimuli

%runs MovBall First



for fi = fileIndex
    fprintf('Starting %s!',files(fi).name)
    tic
    run(fullfile(files(fi).folder, files(fi).name));
    toc
    fprintf('Finished %s!',files(fi).name)
     
end

