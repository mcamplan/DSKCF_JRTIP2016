% runDSKCF.m is a script to test the DS-KCF RGBD tracker described in
%
%[1] S. Hannuna, M. Camplani, J. Hall, M. Mirmehdi, D. Damen, T. Burghardt,
% A.Paiement, L. Tao, DS-KCF: A ~real-time tracker for RGB-D data, Journal
% of Real-Time Image Processing
%
% the script shows how to:
%   -to set DS-KCF parameters
%   -to select folder(s) containing depth and color images 
%   -to run the DS-KCF tracking core
%
%  in this example the RGBD princeton data [2] (validation set) presented in 
%  is used (see the data folder of this package). The data folder structure 
%  has to be as in the example
%  
%  TOPFOLDER
%       VideoSequence1
%               .
%               .
%       VideoSequenceK
%                       depth           (folder Containing Depth Data)
%                       rgb             (folder Containing Color Data)
%                       init.txt        (tracker BB for the first frame)
%                       frames.mat*      (info about frames order)
%                       framesNEW.mat*   (info about frames order)
%               .
%               .
%       VideoSequenceN
%
%  *Please read the DS-KCF paper about depth and color depth stream
%  alignement. framesNEW.mat is the more accurate selection of aligned
%  depth and color frames proposed in our paper for the dataset presented
%  in [2]. Note that to load other datasets or sequences framesNEW.mat could 
%  not be required (i.e. simple sequential images ID)
%
%  About the meaning of the various parameters please read [1] and [3] or
%  the help of WRAPPERDSKCF functions
%
% [1] S. Hannuna, M. Camplani, J. Hall, M. Mirmehdi, D. Damen, T.
% Burghardt, A.Paiement, L. Tao, DS-KCF: A ~real-time tracker for RGB-D
% data, Journal of Real-Time Image Processing
%  
%  [2] S. Song and J. Xiao. Tracking revisited using RGBD camera: Unified benchmark and
%      baselines. In Computer Vision (ICCV), 2013 IEEE International Conference on, pages
%      233–240, 2013.
%
%  [3] J. F. Henriques, R. Caseiro, P. Martins, and J. Batista. High-speed
%  tracking with kernelized correlation filters. Pattern Analysis and
%  Machine Intelligence, IEEE Transactions on, 2015.
%
%  University of Bristol
%  Massimo Camplani and Sion Hannuna
%  
%  massimo.camplani@bristol.ac.uk
%  hannuna@compsci.bristol.ac.uk

clc
clear all
currentFolder=cd();
%%add the DS-KCFresults

dskcfPath{1}='../';
dskcfPath{2}='../functionsDepthSeg';
dskcfPath{3}='../functionsIO';
dskcfPath{4}='../functionsOcclusions';
dskcfPath{5}='../functionsScaleChange';
dskcfPath{6}='../functionsTracking';
dskcfPath{7}='../functionsShape';


for i=1:length(dskcfPath)
    cd(dskcfPath{i});
    tmpPath=cd();
    %addpath(genpath(tmpPath));
    addpath(tmpPath);
    cd(currentFolder);
end

cd(currentFolder)

%insert here the absolute path here you want to save your results or use
%the relative path DS-KCFresults 
rootDestFolder=('../../DS-KCFresults');

mkdir(rootDestFolder);
cd(rootDestFolder);
%take absolute value and create the results folder
rootDestFolder=cd();


cd(currentFolder)

%now select the data folder
rootSourceFolder=('../../data/ValidationSet');
cd(rootSourceFolder);
rootSourceFolder=cd();


%select all the videos in the folder
dirInfo = dir();            
isDir = [dirInfo.isdir];             
listAllVideos = {dirInfo(isDir).name};   
listAllVideos = listAllVideos(3:end);

%If you don't want to precess all the video set this to false
processAllVideos=true;

%eventually select your subset of videos
if(processAllVideos==false)
    %insert video names manually!!!!
    %listVideos{1}='bear_front';
    %listVideos{1}='new_ex_occ4';
    %listVideos{2}='zcup_move_1';
    listVideos{1}='child_no1';    
    %listVideos{1}='face_occ5';    
else
    listVideos=listAllVideos;
end

%Saving info struct saveParameters contains some flags about saving tracker
%results
saveParameters.savingImage=false; %save images with Bounding Box overlapped
saveParameters.savingTrackBool=true; %save tracker output
saveParameters.savingDSKCFParamMat=true;  %save in a matlab file the parameters used
saveParameters.overwriteFolderBool=false; %overwrite results folder or generate new ones
saveParameters.show_visualization=false; %show the tracking results live in a matlab figure
saveParameters.noBitShift=false; %for some dataset the depth data need a bitshift (see wrapperDSKCF())


%% SETTING TRACKER'S PARAMETERS
%  the struct "DSKCFparameters" is built to contains all the parameters it
%  will be created at the end of the section
kernel_type='gaussian';

%change only this flag for feature selection, the rest is automatic!!!!
feature_type = 'hog_concatenate';
kernel.type = kernel_type;

%Different features that can be used
features.rawDepth= false;
features.rawColor=false;
features.rawConcatenate=false;
features.rawLinear=false;
features.hog_color = false;
features.hog_depth = false;
features.hog_concatenate = false;
features.hog_linear = false;


padding = 1.5;  %extra area surrounding the target
lambda = 1e-4;  %regularization
output_sigma_factor = 0.1;  %spatial bandwidth (proportional to target)

%Set the scale Sq in [1]
scales = 0.4:0.1:2.2;


%Note this switch is not necessary, you can eventually 
switch feature_type
    case 'rawDepth',
        interp_factor = 0.075;  %linear interpolation factor for adaptation
        
        kernel.sigma = 0.2;  %gaussian kernel bandwidth
        
        kernel.poly_a = 1;  %polynomial kernel additive term
        kernel.poly_b = 7;  %polynomial kernel exponent
        
        features.rawDepth = true;
        cell_size = 1;
    case 'rawColor',
        interp_factor = 0.075;  %linear interpolation factor for adaptation
        
        kernel.sigma = 0.2;  %gaussian kernel bandwidth
        
        kernel.poly_a = 1;  %polynomial kernel additive term
        kernel.poly_b = 7;  %polynomial kernel exponent
        
        features.rawColor = true;
        cell_size = 1;

      case 'rawConcatenate',
        interp_factor = 0.075;  %linear interpolation factor for adaptation
        
        kernel.sigma = 0.2;  %gaussian kernel bandwidth
        
        kernel.poly_a = 1;  %polynomial kernel additive term
        kernel.poly_b = 7;  %polynomial kernel exponent
        
        features.rawConcatenate = true;
        cell_size = 1;
    case 'rawLinear',
        interp_factor = 0.075;  %linear interpolation factor for adaptation
        
        kernel.sigma = 0.2;  %gaussian kernel bandwidth
        
        kernel.poly_a = 1;  %polynomial kernel additive term
        kernel.poly_b = 7;  %polynomial kernel exponent
        
        features.rawLinear = true;
        cell_size = 1;

    case 'hog_color'
        interp_factor = 0.02;
        
        kernel.sigma = 0.5;
        
        kernel.poly_a = 1;
        kernel.poly_b = 9;
        
        features.hog_color = true;
        features.hog_orientations = 9;
        cell_size = 4;
    case 'hog_depth'
        interp_factor = 0.02;
        
        kernel.sigma = 0.5;
        
        kernel.poly_a = 1;
        kernel.poly_b = 9;
        
        features.hog_depth = true;
        features.hog_orientations = 9;
        cell_size = 4;
    case 'hog_concatenate'
        interp_factor = 0.02;
        
        kernel.sigma = 0.5;
        
        kernel.poly_a = 1;
        kernel.poly_b = 9;
        
        features.hog_concatenate = true;
        features.hog_orientations = 9;
        cell_size = 4;
    case 'hog_linear'
        interp_factor = 0.02;
        
        kernel.sigma = 0.5;
        
        kernel.poly_a = 1;
        kernel.poly_b = 9;
        
        features.hog_linear = true;
        features.hog_orientations = 9;
        cell_size = 4;
                
    otherwise
        error('Unknown feature.')
end

%copy the parameters to the struct
DSKCFparameters.features=features; %feature selection for tracking
DSKCFparameters.kernel=kernel; %kernel size and type
DSKCFparameters.interp_factor=interp_factor; %interpolation factor
DSKCFparameters.cell_size=cell_size; %HOG parameters
DSKCFparameters.padding=padding;
DSKCFparameters.lambda=lambda; 
DSKCFparameters.output_sigma_factor=output_sigma_factor;
DSKCFparameters.scales=scales; % fixed scales

%% PROCESSING LOOP

numVideo=length(listVideos);

% For the Princeton dataset
loadFileOrder='newOrder';
%loadFileOrder='princetonMAT';

%Please read the DS-KCF paper about depth and color depth stream
%  alignement. framesNEW.mat is the more accurate selection of aligned
%  depth and color frames proposed in our paper for the dataset presented
%  in [2]. Note that to load other datasets or sequences framesNEW.mat
%  could not be required (i.e. simple sequential images ID)

%For each selected sequence start to process!!!!!!
for i=1:numVideo
    
    %Generate a tmpResults folder with this name
    %DS-KCF_VIDEONAME_featureSelected_#ID. The folder is generated with a
    %different #ID if run multiple times and
    %saveParameters.overwriteFolderBool has been set to false. See the help
    %of generateFolderResults
    tmpDestFolder=generateFolderResults(rootDestFolder,listVideos{i},feature_type);
    
    %save in the folder input parameters and the current WORKSPACE!!!!!
    save([tmpDestFolder '/inputParam.mat']); 
    
    boolPrinceton=strcmp(loadFileOrder,'princetonMAT');
    boolNewOrder=strcmp(loadFileOrder,'newOrder');
    
    %Load Image information as path etc. Read also the inititial target
    %position
    
    if(boolNewOrder)
        [img_files, depth_files, pos, target_sz, ground_truth, video_path, depth_path] = ...
            load_video_info_depthFROMMAT(rootSourceFolder, listVideos{i},1);
    else
        [img_files, depth_files, pos, target_sz, ground_truth, video_path, depth_path] = ...
            load_video_info_depthFROMMAT(rootSourceFolder, listVideos{i},0);
    end
    
    %call tracker wrapper function with all the relevant parameters
    [dsKCFoutputSr,dsKCFoutputSq,dsKCFsegmentationOut, avTime,totalTime,timeMatrix] = ...
        wrapperDSKCF(video_path, depth_path,img_files, depth_files, pos, ...
        target_sz, DSKCFparameters,saveParameters.show_visualization,...
        saveParameters.savingImage,tmpDestFolder,saveParameters.noBitShift );
   
    %save tracking results and processing time
    save([tmpDestFolder '/procTime.mat'], 'totalTime');
    save([tmpDestFolder '/procTime.txt'], 'totalTime','-ascii');
    
    save([tmpDestFolder '/procTimeFrames.mat'], 'avTime');
    save([tmpDestFolder '/procTimeFrames.txt'], 'avTime','-ascii');

    save([tmpDestFolder '/modulesTimeFrames.mat'], 'timeMatrix');
    save([tmpDestFolder '/modulesTimeFrames.txt'], 'timeMatrix','-ascii');

    
    %Note we are saving the results as y x height width to be consistent
    %with the notation presented in the princeton RGB-D [2]

    %Results using Sr in [1] use this for your comparison
    trackRes=[dsKCFoutputSr];
    save([tmpDestFolder '/' listVideos{i} '.mat'], 'trackRes');
    save([tmpDestFolder '/' listVideos{i} '.txt'], 'trackRes','-ascii');

    %Results using Sq in [1]
    trackResTargetSize=[dsKCFoutputSq];
    save([tmpDestFolder '/' listVideos{i} 'TargetSize.mat'], 'trackResTargetSize');
    save([tmpDestFolder '/' listVideos{i} 'TargetSize.txt'], 'trackResTargetSize','-ascii');
    
    %Results using segmentation
    save([tmpDestFolder '/' listVideos{i} 'dsKCFsegmentationOut.mat'], 'dsKCFsegmentationOut');
    save([tmpDestFolder '/' listVideos{i} 'dsKCFsegmentationOut.txt'], 'dsKCFsegmentationOut','-ascii');
end