
%Set up python env to be able to read pickle
pyenv('Version', 'C:\Users\MarkS9\anaconda3\envs\py9\python.exe');
%%
csvPaths = 'C:\Users\MarkS9\Desktop\HeadFixedEyeTrack2-Simon-2024-10-02\EyeTracking3-Simon-2024-10-03\videos';
% Import the Python pickle module
%pickle = py.importlib.import_module('pickle');

% Use dir to find files containing 'full.pickle' in their names
fileList = dir(fullfile(csvPaths, '*.csv*'));

% Extract the file names into a cell array
fileNames = {fileList.name}';

% %% Read pickle
% % Open the pickle file
% fid = py.open([picklePaths filesep fileNames{1}]);
% % Load the data from the pickle file
% data = pickle.load(fid);
% % Close the file
% fid.close();
% % Convert the data to a MATLAB variable
% matlab_data = py2mat(data);
% %picklePaths = '/C:/Users/MarkS9/Desktop/HeadFixedEyeTrack2-Simon-2024-10-02/EyeTracking3-Simon-2024-10-03/videos';
% % Use dir to find files containing 'full.pickle' in their names
% fileList = dir(fullfile(picklePaths, '*snapshot_200.h5*'));
% %Extract the file names into a cell array
% fileNames = {fileList.name}';
% data = h5read(fileNames{1},picklePaths);
% cd(picklePaths)
% %%
% 
% % I=h5info([picklePaths filesep fileNames{1},'/'])
% % data = h5read([picklePaths filesep fileNames{1}],'/g4/abounds');
% 
% 
% %%
% 
% dataP = py.pickle.load(py.open(fullfile(picklePaths,fileNames{1}),'rb'));
% 
% dataPm = sruct(dataP);


% Read the entire CSV file into a table
Data = readtable([csvPaths filesep fileNames{1}]);



