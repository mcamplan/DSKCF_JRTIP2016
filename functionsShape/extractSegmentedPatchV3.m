function [ reshapedSegmentedMask, insidePatchIndexes, finalSegmentedMask] = ...
    extractSegmentedPatchV3( segmentedMask,outOfBoundSize,centralPoint, shapeDSKCF_struct,maxSize )
% EXTRACTSEGMENTEDPATCHV3.m function manage the alignment of segmented
% patches to accumulate object binary masks as in [1]
%
%
%   EXTRACTSEGMENTEDPATCHV3 is used to accumulate the segmented target
%   masks and align them correctly in case of target close to the borders
%   for example or change in scale while those silhouette are accumulated
%
%   INPUT:
%   - segmentedMask   segmented binary mask of the current frame
%   - outOfBoundSize  size of the patch to be assigned if the vector
%   containing the segmented masks is empty
%   - centralPoint center of the current tracking target
%   -shapeDSKCF_struct data structure containing shape information (see INITDSKCFSHAPE) 
%   - maxSize image size
%   OUTPUT
%   -reshapedSegmentedMask segmented mask with adjusted size according to
%   the accumulated mask vector
%   -insidePatchIndexes index of the mask inside the patch (needed to be
%   re-adjusted in case of object close to borders of change of scale)see
%   also ADDSEGMENTATIONRESULTS
%   -finalSegmentedMask a copy of input segmentedMask
%
%  See also INITDSKCFSHAPE, ADDSEGMENTATIONRESULTS
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

stPointX=1;
stPointY=1;
endPointY=size(segmentedMask,1);
endPointX=size(segmentedMask,2);

stPointXMask=1;
stPointYMask=1;
endPointYMask=size(segmentedMask,1);
endPointXMask=size(segmentedMask,2);

%clip the groundflor
segmentedMask(round(endPointYMask*0.85):end,round(endPointXMask*0.60):end)=0;
segmentedMask(round(endPointYMask*0.85):end,1:round(endPointXMask*0.4))=0;

if(isempty(shapeDSKCF_struct.cumulativeMask))
    cumulativeSize=outOfBoundSize;
else
    cumulativeSize=size(shapeDSKCF_struct.cumulativeMask);
end

if(size(segmentedMask)==cumulativeSize)
    reshapedSegmentedMask=segmentedMask;
else
    extremaBB(1:2)=centralPoint-cumulativeSize(2:-1:1)/2;
    extremaBB(3:4)=centralPoint+cumulativeSize(2:-1:1)/2;
    extremaBB=floor(extremaBB(:));
    reshapedSegmentedMask=repmat(0,cumulativeSize);
    
    %CASE1: object on the left
    if(extremaBB(1)<0)%problem with x placement
        diffSize=min(size(segmentedMask,2)-cumulativeSize(2),0);
        %%dif size is negative or zero
        stPointX=max(1,abs(extremaBB(1))+diffSize);
        minSize=min(size(segmentedMask,2),cumulativeSize(2));
        endPointX=stPointX+minSize-1;
        stPointXMask=1;
        endPointXMask=stPointXMask+minSize-1;
    end
    
    %CASE1: object on the right
    if(extremaBB(3)>maxSize(1))
        minSize=min(size(segmentedMask,2),cumulativeSize(2));
        offset=extremaBB(3)-maxSize(1);
        endPointX=cumulativeSize(2)-offset;
        stPointX=max((endPointX-minSize+1),1);
        endPointXMask=minSize;
        numberOfFinalElement=endPointX-stPointX;
        stPointXMask=endPointXMask-numberOfFinalElement;
    end
    
    %CASE3: object on the top
    if(extremaBB(2)<0)%problem with x placement
        
        diffSize=min(size(segmentedMask,1)-cumulativeSize(1),0);
        stPointY=max(1,abs(extremaBB(2))+diffSize);
        minSize=min(size(segmentedMask,1),cumulativeSize(1));
        endPointY=stPointY+minSize-1;
        stPointYMask=1;
        endPointYMask=stPointYMask+minSize-1;
    end
    
    %CASE1: object on the bottom
    if(extremaBB(4)>maxSize(2))
        minSize=min(size(segmentedMask,1),cumulativeSize(1));
        offset=extremaBB(4)-maxSize(2);
        endPointY=cumulativeSize(1)-offset;
        stPointY=max((endPointY-minSize+1),1);
        numberOfFinalElement=endPointY-stPointY;
        stPointYMask=endPointYMask-numberOfFinalElement;
    end
    
    reshapedSegmentedMask(stPointY:endPointY,stPointX:endPointX)=segmentedMask(stPointYMask:endPointYMask,stPointXMask:endPointXMask);
    
end

selectedSize=[length(stPointY:endPointY),length(stPointX:endPointX)];

minSIZE=[min(size(segmentedMask,1),cumulativeSize(1)), min(size(segmentedMask,2),cumulativeSize(2))];


if(selectedSize(2)<minSIZE(2))%cumulativeSize(2))
    diffLeft=stPointX-1;
    diffRight=minSIZE(2)-endPointX;%cumulativeSize(2)-endPointX;
    endPointX=endPointX+diffRight;
    stPointX=stPointX-diffLeft;
end

if(selectedSize(1)<minSIZE(1) )%cumulativeSize(1))
    diffUp=stPointY-1;
    diffDown=minSIZE(1)-endPointY;%cumulativeSize(1)-endPointY;
    endPointY=endPointY+diffDown;
    stPointY=stPointY-diffUp;
end

selectedSizeFinal=size(reshapedSegmentedMask);

if(selectedSizeFinal(2)>cumulativeSize(2))%)
    diff=selectedSizeFinal(2)-minSIZE(2);
    diffLIMIT=selectedSizeFinal(2)-cumulativeSize(2);
    endPointXLIMIT=endPointX-diffLIMIT;
    endPointX=endPointX-diff;
    stPointX=1;
    reshapedSegmentedMask=reshapedSegmentedMask(:,stPointX:endPointXLIMIT);
end

if(selectedSizeFinal(1)>cumulativeSize(1))%cumulativeSize(2))
    diff=selectedSizeFinal(1)-minSIZE(1);
    diffLIMIT=selectedSizeFinal(1)-cumulativeSize(1);
    endPointYLIMIT=endPointY-diffLIMIT;
    endPointY=endPointY-diff;
    stPointY=1;
    reshapedSegmentedMask=reshapedSegmentedMask(stPointY:endPointYLIMIT,:);
end


insidePatchIndexes=[stPointX,stPointY,endPointX,endPointY];


finalSegmentedMask=segmentedMask;
end
