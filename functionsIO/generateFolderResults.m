% generateFolderResults.m is a function to automatically generate the results folder where tracker's output files are stored
%
%  GENERATEFOLDERRESULTS creates the results folder for the DS-KCF tracker
%  inside the top folder specified in the string rootDestFolder. The
%  generated folder name is composed by using two other strings: videoName
%  that contains the sequence name that is going to be analyzed and
%  feature_type string that contains the name of the feature selected for
%  the DS-KCF tracker. In the case that a folder with the same name already
%  exists in rootDestFolder, an incremental counter is used to avoid
%  overwriting the old folder
%
%   INPUT:
%  -rootDestFolder name of the top folder where results will be saved
%  -videoName name of the sequences processed by the DS-KCF tracker
%  -feature_type name of the feauture used by the DS-KCF tracker
%
%   OUTPUT
%  -tmpDestFolder string containing absolute path of the results folder
%
%  Examples:
%    >> rootDestFolder='C:\myExistingResultsFolder';
%    >> videoName='videoToProcess';
%    >> feature_type='hog_depth';
%    >> tmpDestFolder=generateFolderResults(rootDestFolder,videoName,feature_type)
%       tmpDestFolder =
%
%       C:\myExistingResultsFolder\videoToProcess_hog_depth
%
%  University of Bristol
%  Massimo Camplani and Sion Hannuna
%
%  massimo.camplani@bristol.ac.uk
%  hannuna@compsci.bristol.ac.uk

function tmpDestFolder=generateFolderResults(rootDestFolder,videoName,feature_type)

tmpDestFolderAbsolute=[rootDestFolder '/' videoName '_' feature_type];
tmpDestFolder=[videoName '/' feature_type];

existingDir=exist(tmpDestFolderAbsolute,'dir');

if(existingDir==false)
    
    tmpDestFolder=tmpDestFolderAbsolute;
else
    %take the number of existing folder...
    numDir=length(dir([tmpDestFolderAbsolute '*']));
    numDir=numDir+1;
    tmpDestFolderAbsolute=[rootDestFolder '/' videoName '_' feature_type num2str(numDir)];
    tmpDestFolder=tmpDestFolderAbsolute;
end

mkdir(tmpDestFolder);