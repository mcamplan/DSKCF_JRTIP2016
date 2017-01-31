function trackerDSKCF_struct=initDSKCFtracker()
% INITDSKCFTRACKER.m initializes the data structure for DS-KCF tracker [1]
% 
%   INITDSKCFTRACKER function initializes the data structure of the
%   DS-KCF tracker. In particular
%
%   INPUT: none
%
%   OUTPUT
%   -trackerDSKCF_struct data structure that contains DS-KCF tracker data
%  structure
%   
%   + currentTarget.posX column in the image plane
%   + currentTarget.posY row in the image plane
%   + currentTarget.h height of the target
%   + currentTarget.w width of the target
%   + currentTarget.bb bounding box of the target in the format 
%                     [topLeftX, topLeftY, bottomRightX, bottomRightY]
%   + currentTarget.meanDepthObj mean depth of the tracker object
%   + currentTarget.stdDepthObj depth's standard deviation of the tracker object
%   + currentTarget.LabelRegions cluster labels of the segmented target region
%   + currentTarget.regionIndex= label of the object cluster
%   + currentTarget.Centers depth centers of the clusters
%   + currentTarget.LUT=[] LUT
%   + currentTarget.occBB=[0 0 0 0]; occluding bounding box in the format
%   [topLeftX, topLeftY, bottomRightX, bottomRightY]
%   + currentTarget.totalOcc=0;  total occlusion flag
%   + currentTarget.underOcclusion=0;  under occlusion flag
%   +currentTarget.conf maximum response of the DSKCF for the current frame
%
%   models in the frequency domain for the KCFbased tracking by using color
%   and depth features (see [1] for mor details)
%   +model_alphaf = []; 
%   +model_alphaDf = [];
%   +model_xf = [];
%   +model_xDf = [];
%
%   +previousTarget contains same information of currentTarget, but they
%   it is relative to the target tracked in the previous frame.
%
%  See also WRAPPERDSKCF, INITDSKCFTRACKER_occluder
%
% [1] S. Hannuna, M. Camplani, J. Hall, M. Mirmehdi, D. Damen, T.
% Burghardt, A.Paiement, L. Tao, DS-KCF: A ~real-time tracker for RGB-D
% data, Journal of Real-Time Image Processing
%
%
%  University of Bristol 
%  Massimo Camplani and Sion Hannuna
%  
%  massimo.camplani@bristol.ac.uk 
%  hannuna@compsci.bristol.ac.uk

trackerDSKCF_struct=[];

% current target position and bounding box
trackerDSKCF_struct.currentTarget.posX=0;%column in the image plane
trackerDSKCF_struct.currentTarget.posY=0;%row in the image plane
trackerDSKCF_struct.currentTarget.h=0;%height of the target
trackerDSKCF_struct.currentTarget.w=0;%width in the image planeof the target
trackerDSKCF_struct.currentTarget.bb=[0 0 0 0]; % in the format [topLeftX, topLeftY, bottomRightX, bottomRightY]
trackerDSKCF_struct.currentTarget.conf=0;
%current target depth distribution info
trackerDSKCF_struct.currentTarget.meanDepthObj=0;% mean depth of the tracker object
trackerDSKCF_struct.currentTarget.stdDepthObj=0;% depth's standard deviation of the tracker object
trackerDSKCF_struct.currentTarget.LabelRegions=[];%cluster labels of the segmented target region
trackerDSKCF_struct.currentTarget.regionIndex=0;%label of the object cluster
trackerDSKCF_struct.currentTarget.Centers=[];%depth centers of the clusters
trackerDSKCF_struct.currentTarget.LUT=[];%LUT
trackerDSKCF_struct.currentTarget.segmentedBB=[];%bounding box of the corresponding sgmented region
%current target depth occluding info
trackerDSKCF_struct.currentTarget.occBB=[0 0 0 0]; % in the format [topLeftX, topLeftY, bottomRightX, bottomRightY]
trackerDSKCF_struct.currentTarget.totalOcc=0; % total occlusion flag
trackerDSKCF_struct.currentTarget.underOcclusion=0; % under occlusion flag
%target model alpha and X, see [1] for more details
trackerDSKCF_struct.model_alphaf = []; 
trackerDSKCF_struct.model_alphaDf = [];
trackerDSKCF_struct.model_xf = [];
trackerDSKCF_struct.model_xDf = [];


%previous target entries
trackerDSKCF_struct.previousTarget=trackerDSKCF_struct.currentTarget;


