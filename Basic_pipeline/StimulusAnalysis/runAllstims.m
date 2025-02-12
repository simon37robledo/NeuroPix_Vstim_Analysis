%%Run all stimuli scripts

files = dir('D:\Mark_S13\Documents\GitHub\NeuroPix_Vstim_Analysis\Basic_pipeline\StimulusAnalysis\*.m'); % Get all .m files
for fi = 1:length(files)
    fprintf('Starting %s!',files(fi).name)
    tic
    run(fullfile(files(fi).folder, files(fi).name));
    toc
    fprintf('Finished %s!',files(fi).name)
end
