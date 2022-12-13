

overallPath = '/Users/matttang/IDrive-Sync/Matt/Work/Current work/Rodent V1 and prediction/' % where the scripts are stored
dataPath='/Users/matttang/Downloads/Saved' % where the data is stored

outputPath=[overallPath 'Output'];
addpath(genpath([ overallPath 'FinalScripts/Functions']));
figPath =[overallPath 'Figures'];
cd(dataPath);
fileNames = dir('*-SaveNoEpoch.mat');