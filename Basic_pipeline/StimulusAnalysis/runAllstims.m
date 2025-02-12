%%Run all stimuli scripts
SDGrecordingsA = [8:14,40:43,49:54]; %anesthetized 
files = dir('D:\Mark_S13\Documents\GitHub\NeuroPix_Vstim_Analysis\Basic_pipeline\StimulusAnalysis\*.m'); 
filenames = {files.name};
fileIndex = find(~contains(filenames,'runAllstims'));% Avoid infinite loop
for fi = fileIndex
    fprintf('Starting %s!',files(fi).name)
    tic
    run(fullfile(files(fi).folder, files(fi).name));
    toc
    fprintf('Finished %s!',files(fi).name)
     
end

