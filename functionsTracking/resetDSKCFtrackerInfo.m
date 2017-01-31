function trackerDSKCF_struct=resetDSKCFtrackerInfo(trackerDSKCF_struct)
% RESETDSKCFTRACKERINFO.m re-initializes the data structure for DS-KCF tracker before processing a new frame[1]
% 
%   RESETDSKCFTRACKERINFO re-initializes the data structure for DS-KCF
%   tracker before processing a new frame[1]
%
%   INPUT: 
%   -trackerDSKCF_structIN tracker data structure
%   OUTPUT
%  -trackerDSKCF_struct data structure with re-set tracking info
%   
%  [1] S. Hannuna, M. Camplani, J. Hall, M. Mirmehdi, D. Damen, T.
%  Burghardt, A. Paiement, L. Tao, DS-KCF: A real-time tracker for RGB-D
%  data, Journal of Real-Time Image Processing
% See also SINGLEFRAMEDSKCF
%
%  University of Bristol 
%  Massimo Camplani and Sion Hannuna
%  
%  massimo.camplani@bristol.ac.uk 
%  hannuna@compsci.bristol.ac.uk

%trackerDSKCF_struct=trackerDSKCF_structIN;
trackerDSKCF_struct.currentTarget.occBB=[0 0 0 0]; % in the format [topLeftX, topLeftY, bottomRightX, bottomRightY]
trackerDSKCF_struct.currentTarget.totalOcc=0; % total occlusion flag
trackerDSKCF_struct.currentTarget.underOcclusion=0; % under occlusion flag
trackerDSKCF_struct.currentTarget.conf=0;

