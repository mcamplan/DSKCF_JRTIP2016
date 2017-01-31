function [trackerDSKCF_structOccluder,newPos]=singleFrameDSKCF_occluder(firstFrame,...
    im,depth,trackerDSKCF_structOccluder,DSKCFparameters)
% SINGLEFRAMEDSKCF_OCCLUDER.m functions for tracking occluding object
%
%   SINGLEFRAMEDSKCF_OCCLUDER is the function for tracking the occluding
%   object in the DS-KCF tracker framework (for more details see [1]). This
%   function is based on several data structures where input and output
%   data is stored.Please note that data structures are supposed to be
%   initialized as in wrapperDSKCF and runDSKCF.m test script.
%
%   INPUT:
%   - firstFrame   boolean flag that marks the first frame of the
%   sequence
%   - im    color data
%   - depth    depth data (8bit image)
%   - trackerDSKCF_structOccluder  DS-KCF tracker data structure (see WRAPPERDSKCF,
%   INITDSKCFTRACKER)
%   - DSKCFparameters, parameters structures
%
%
%   OUTPUT
%   -newPos updated position of the DS-KCF tracker while tracking the
%   occluding object newPos=[y x] where x is the column and y is the row
%   index
%   -trackerDSKCF_structOccluder updated data structure for the trackers
%
%  See also MAXRESPONSEDSKCF , SINGLEFRAMEDSKCF, GET_SUBWINDOW,
%  MODELUPDATEDSKCF, FROMCENTRALPOINTTOBB
%
%
%  [1] S. Hannuna, M. Camplani, J. Hall, M. Mirmehdi, D. Damen, T.
%  Burghardt, A. Paiement, L. Tao, DS-KCF: A real-time tracker for RGB-D
%  data, Journal of Real-Time Image Processing
%
%
%  University of Bristol
%  Massimo Camplani and Sion Hannuna
%
%  massimo.camplani@bristol.ac.uk
%  hannuna@compsci.bristol.ac.uk


%insert in pos the previous target position...
pos=[trackerDSKCF_structOccluder.previousTarget.posY,trackerDSKCF_structOccluder.previousTarget.posX];
newPos=pos;
if(firstFrame==false)
    
    %obtain a subwindow for training at newly estimated target position
    patch = get_subwindow(im, pos, trackerDSKCF_structOccluder.window_sz);
    patch_depth = get_subwindow(depth, pos, trackerDSKCF_structOccluder.window_sz);
    nRows=size(im,1);
    nCols=size(im,2);
    %calculate response of the DS-KCF tracker
    [response, maxResponse,newPos]=maxResponseDSKCF(...
        patch,patch_depth, DSKCFparameters.features,DSKCFparameters.kernel,...
        pos,DSKCFparameters.cell_size, trackerDSKCF_structOccluder.cos_window,...
        trackerDSKCF_structOccluder.model_xf,trackerDSKCF_structOccluder.model_alphaf,...
        trackerDSKCF_structOccluder.model_xDf,trackerDSKCF_structOccluder.model_alphaDf,...
        nRows,nCols);
        
    %update tracker struct, new position etc
    trackerDSKCF_structOccluder.currentTarget.posX=newPos(2);
    trackerDSKCF_structOccluder.currentTarget.posY=newPos(1);
    
    trackerDSKCF_structOccluder.currentTarget.bb=fromCentralPointToBB...
        (trackerDSKCF_structOccluder.currentTarget.posX,trackerDSKCF_structOccluder.currentTarget.posY,...
        trackerDSKCF_structOccluder.currentTarget.w,trackerDSKCF_structOccluder.currentTarget.h,...
        nCols,nRows);
    
    trackerDSKCF_structOccluder.currentTarget.conf=maxResponse;
    
end

%obtain a subwindow for training at newly estimated target position
patch = get_subwindow(im, newPos, trackerDSKCF_structOccluder.window_sz);
patch_depth = get_subwindow(depth, newPos, trackerDSKCF_structOccluder.window_sz);

%update occluder model....
[trackerDSKCF_structOccluder.model_alphaf, trackerDSKCF_structOccluder.model_alphaDf, ...
    trackerDSKCF_structOccluder.model_xf, trackerDSKCF_structOccluder.model_xDf]=...
    modelUpdateDSKCF(firstFrame,patch,patch_depth,DSKCFparameters.features,...
    DSKCFparameters.cell_size,trackerDSKCF_structOccluder.cos_window,...
    DSKCFparameters.kernel,trackerDSKCF_structOccluder.yf,...
    DSKCFparameters.lambda,trackerDSKCF_structOccluder.model_alphaf, ...
    trackerDSKCF_structOccluder.model_alphaDf,trackerDSKCF_structOccluder.model_xf,...
    trackerDSKCF_structOccluder.model_xDf,0,DSKCFparameters.interp_factor);

