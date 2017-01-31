function [ dskcfShapeStruct ] = initDSKCFshape( slidingWinSize, scaleFactor,dskcfShapeStructOLD )
% INITDSKCFSHAPE.m initializes the shape data structure for DS-KCF tracker
% presented in [1]
%
%   INITDSKCFSHAPE function initializes the shape data structure of the
%   DS-KCF tracker. In particular, some matrices are precomputed at the
%   different scales and other flags are intialized
%
%   INPUT:
%  -slidingWinSize size of the sliding window used to accumulate silhouette
%  masks as in [1]
%  -current scaleFactor of the DSKCF, when zero the dskcfShapeStruct is
%  initialized with empty vectors
%  -pos initial target position of the tracked object
%   OUTPUT
%  -dskcfShapeStructOLD data structure that contains shape parameters to be
%  updated. It contains
%
%  + slidingWindowSize size of the sliding windows to accumulated segmented
%  masks
%  + lastSegmentedBB bounding box containing the segmented object
%  + cumulativeMask mask obtained accumulating masks in the sliding window
%  as in [1]
%  +maskArray array containing the segmented masks
%  +segmentW width of the segmented object
%  +segmentH height of the segmented object
%  +growingStatus flag to mark  if the object is increasing its size while
%  changing shape
%
%  [1] S. Hannuna, M. Camplani, J. Hall, M. Mirmehdi, D. Damen, T.
%  Burghardt, A.Paiement, L. Tao, DS-KCF: A real-time tracker for RGB-D
%  data, Journal of Real-Time Image Processing
%
%
%  University of Bristol
%  Massimo Camplani and Sion Hannuna
%
%  massimo.camplani@bristol.ac.uk
%  hannuna@compsci.bristol.ac.uk

dskcfShapeStruct.slidingWindowSize=slidingWinSize;
dskcfShapeStruct.lastSegmentedBB=[];

if(scaleFactor>0)
    
    dskcfShapeStruct.cumulativeBB=dskcfShapeStructOLD.cumulativeBB;
    dskcfShapeStruct.growingStatus=dskcfShapeStructOLD.growingStatus;
    
    %zero pad previous segmentations
    newW=round(scaleFactor*dskcfShapeStructOLD.segmentW);
    newH=round(scaleFactor*dskcfShapeStructOLD.segmentH);
    
    segmentHIncrement=round((newH-dskcfShapeStructOLD.segmentH)/2);
    segmentWIncrement=round((newW-dskcfShapeStructOLD.segmentW)/2);
    
    dskcfShapeStruct.segmentH=newH;
    dskcfShapeStruct.segmentW=newW;
    
    segmentHIncrement=segmentHIncrement*(segmentHIncrement>0);
    segmentWIncrement=segmentWIncrement*(segmentWIncrement>0);
    
    
    dskcfShapeStruct.cumulativeMask=padarray(dskcfShapeStructOLD.cumulativeMask,[segmentHIncrement segmentWIncrement]);
    %you need only from the second one...
    
    for i=1:size(dskcfShapeStructOLD.maskArray,3)
        tmpMaskArray(:,:,i)=...
            padarray(dskcfShapeStructOLD.maskArray(:,:,i),[segmentHIncrement segmentWIncrement]);
    end
    dskcfShapeStruct.maskArray=tmpMaskArray;
    
else
    dskcfShapeStruct.cumulativeMask=[];
    dskcfShapeStruct.maskArray=[];
    dskcfShapeStruct.cumulativeBB=[];
    dskcfShapeStruct.segmentW=0;
    dskcfShapeStruct.segmentH=0;
    dskcfShapeStruct.growingStatus=false;
end

end

