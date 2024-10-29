
%Set up python env to be able to read pickle
pyenv('Version', 'C:\Users\MarkS9\anaconda3\envs\py9\python.exe');

picklePaths = '/C:\Users\MarkS9\Desktop\HeadFixedEyeTrack2-Simon-2024-10-02\EyeTracking3-Simon-2024-10-03\videos';

picklePaths = '/C:/Users/MarkS9/Desktop/HeadFixedEyeTrack2-Simon-2024-10-02/EyeTracking3-Simon-2024-10-03/videos';

% Use dir to find files containing 'full.pickle' in their names
fileList = dir(fullfile(picklePaths, '*snapshot_200.h5*'));

% Extract the file names into a cell array
fileNames = {fileList.name}';

data = h5read(fileNames{1},picklePaths);

cd(picklePaths)



%%

dataP = py.pickle.load(py.open(fullfile(picklePaths,fileNames{1}),'rb'));

dataPm = sruct(dataP);