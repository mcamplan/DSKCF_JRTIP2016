function [estimatedShapeBB,shapeDSKCF_struct,changeOfShapeFlag,newOutput]=...
    regionModificationCheck(sizeOfSegmenter,sizeOfTarget,accumulatedSEGBool...
    ,noDataPercent,minSizeOK,estimatedShapeBB,shapeDSKCF_struct,imageSize...
    ,trackerDSKCF_struct)
% REGIONMODIFICATIONCHECK.m function to identify change of shape as
% presented in [1]
%
%
%   REGIONMODIFICATIONCHECK is used to estimate significant change of shape
%   of the segmented target silhouette, if detected the new shape is used
%   as output of the DS-KCF tracker as presented in [1]
%
%   INPUT:
%   - sizeOfSegmenter   size of the segmented object
%   - sizeOfTarget  size of the tracked target
%   - accumulatedSEGBool binary flag to mark if the vector containing the
%   segmented masks (see [1]) is full, so the segmentation is considered
%   reliable
%   -noDataPercent binary flag to identify if the percentage of missing
%   depth data is greater than a threshold (set in SINGLEFRAMEDSKCF). In
%   case that the flag has value 1 segmentation is not considered reliable
%   -minSizeOK binary flag to identify if the size of the object is greater
%   than a threshold (set in SINGLEFRAMEDSKCF). In case that the flag has
%   value 1 segmentation is not considered reliable
%   -estimatedShapeBB boundinbox of the segmented object obtained with the
%   accumulation strategy presented in [1]
%   -shapeDSKCF_struct data structure containing shape information (see
%   INITDSKCFSHAPE)
%   - imSize image size
%   - trackerDSKCF_struct  DS-KCF tracker data structure (see WRAPPERDSKCF,
%   INITDSKCFTRACKER)
%
%   OUTPUT -estimatedShapeBB modified BB of the segmented target
%   -shapeDSKCF_struct modified data structure containing shape information
%   (see INITDSKCFSHAPE) -changeOfShapeFlag,newOutput flags to mark a
%   detected change of shape according to the method presented in [1]
%
%  See also INITDSKCFSHAPE, EXTRACTSEGMENTEDPATCHV3
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

%simple case, the segmentar has not grown yet, check for  smaller object or
%a bigger one
newOutput=false;
changeOfShapeFlag=false;
if(accumulatedSEGBool && noDataPercent)
    if(shapeDSKCF_struct.growingStatus==false)
        
        if(sizeOfSegmenter<sizeOfTarget*0.9  && minSizeOK )
            estimatedShapeBB=enlargeBB(estimatedShapeBB ,-0.05,imageSize);
            shapeDSKCF_struct.growingStatus=false;
            newOutput=true;
        end
        
        if(sizeOfSegmenter>sizeOfTarget*1.09 )
            %trackerDSKCF_struct.currentTarget.segmentedBB=tmpBBforSegCumulative(:)';
            changeOfShapeFlag=true;
            newOutput=true;
        end
    else
        
        [centerX,centerY,width,height]=fromBBtoCentralPoint(estimatedShapeBB);
        estimatedShapeSize=[height,width];
        
        currentSizeSegmenter=[shapeDSKCF_struct.segmentH,shapeDSKCF_struct.segmentW];
        currentSizeMASK=[size(shapeDSKCF_struct.maskArray(:,:,1),1),size(shapeDSKCF_struct.maskArray(:,:,1),2)];
        
        diffSize=estimatedShapeSize-currentSizeSegmenter;
        
        growingFlag=(diffSize(1)>0 || diffSize(2)>0);
        shrinkingFlag=(diffSize(1)<=0 && diffSize(2)<=0);
        %if continues to grow or decrease the size ....then reshape....
        
        if(growingFlag)
            segmentHIncrement=0;
            segmentWIncrement=0;
            newOutput=true;
            if(diffSize(1)>0.07*currentSizeSegmenter(1))
                segmentHIncrement=round(0.05*currentSizeMASK(1));
                shapeDSKCF_struct.growingStatus=true;
                %newOutput=true;
            end
            
            if(diffSize(2)>0.05*currentSizeSegmenter(2))
                segmentWIncrement=round(0.05*currentSizeMASK(2));
                shapeDSKCF_struct.growingStatus=true;
                %newOutput=true;
            end
            
            shapeDSKCF_struct.growingStatus=true;
            if(segmentHIncrement>0 || segmentWIncrement>0)
                shapeDSKCF_struct.cumulativeMask=padarray(shapeDSKCF_struct.cumulativeMask,[segmentHIncrement segmentWIncrement]);
                for i=1:size(shapeDSKCF_struct.maskArray,3)
                    tmpMaskArray(:,:,i)=...
                        padarray(shapeDSKCF_struct.maskArray(:,:,i),[segmentHIncrement segmentWIncrement]);
                end
                shapeDSKCF_struct.maskArray=tmpMaskArray;
                
                shapeDSKCF_struct.segmentW=max([estimatedShapeSize(2),trackerDSKCF_struct.currentTarget.w]);
                shapeDSKCF_struct.segmentH=max([estimatedShapeSize(1),trackerDSKCF_struct.currentTarget.h]);
            end
            
        end
        %else if shrink a lot....back to the normal status
        if(shrinkingFlag)
            
            %find the cropping area...
            %is inside the target BB?
            if(estimatedShapeSize(1)<trackerDSKCF_struct.currentTarget.h ...
                    && estimatedShapeSize(2)<trackerDSKCF_struct.currentTarget.w)
                
                bbIn=fromCentralPointToBB(round(currentSizeMASK(2)/2),round(currentSizeMASK(1)/2),...
                    estimatedShapeSize(2),estimatedShapeSize(1),...
                    currentSizeMASK(2),currentSizeMASK(1));
                
                
                newOutput=false;%decide later what you want to show
                shapeDSKCF_struct.growingStatus=false;
                
                shapeDSKCF_struct.segmentW=trackerDSKCF_struct.currentTarget.w;
                shapeDSKCF_struct.segmentH=trackerDSKCF_struct.currentTarget.h;
                
                
                %reduce the search region
                %grow the mask properly...
                
                tmpCumulativeMask=roiFromBB(shapeDSKCF_struct.cumulativeMask,bbIn);
                segmentHIncrement=round(((1.1*trackerDSKCF_struct.currentTarget.h)-size(tmpCumulativeMask,1))/2);
                
                segmentWIncrement=round(((1.1*trackerDSKCF_struct.currentTarget.w)-size(tmpCumulativeMask,2))/2);
                %be sure >0
                segmentHIncrement=segmentHIncrement*(segmentHIncrement>0);
                segmentWIncrement=segmentWIncrement*(segmentWIncrement>0);
                
                shapeDSKCF_struct.cumulativeMask=padarray(tmpCumulativeMask,[segmentHIncrement segmentWIncrement]);
                for i=1:size(shapeDSKCF_struct.maskArray,3)
                    tmpMaskArray(:,:,i)=...
                        roiFromBB(shapeDSKCF_struct.maskArray(:,:,i),bbIn);
                end
                shapeDSKCF_struct.maskArray=padarray(tmpMaskArray,[segmentHIncrement segmentWIncrement]);
                
                
                
            else
                
                bbIn=fromCentralPointToBB(round(currentSizeMASK(2)/2),round(currentSizeMASK(1)/2),...
                    estimatedShapeSize(2),estimatedShapeSize(1),...
                    currentSizeMASK(2),currentSizeMASK(1));
                
                newOutput=true;%show this one!!!!!
                shapeDSKCF_struct.growingStatus=true;
                
                shapeDSKCF_struct.segmentW=max([estimatedShapeSize(2),trackerDSKCF_struct.currentTarget.w]);
                shapeDSKCF_struct.segmentH=max([estimatedShapeSize(1),trackerDSKCF_struct.currentTarget.h]);
                
                
                %reduce the search region
                %shrink the masks, but eventually you need to repad the
                %area as the minimun target area
                tmpCumulativeMask=roiFromBB(shapeDSKCF_struct.cumulativeMask,bbIn);
                segmentHIncrement=round(((trackerDSKCF_struct.currentTarget.h)-estimatedShapeSize(1))/2);
                %segmentWIncrement=round(0.05*trackerDSKCF_struct.currentTarget.w);
                segmentWIncrement=round(((trackerDSKCF_struct.currentTarget.w)-estimatedShapeSize(2))/2);
                %be sure >0
                segmentHIncrement=segmentHIncrement*(segmentHIncrement>0);
                segmentWIncrement=segmentWIncrement*(segmentWIncrement>0);
                
                shapeDSKCF_struct.cumulativeMask=padarray(tmpCumulativeMask,[segmentHIncrement segmentWIncrement]);
                for i=1:size(shapeDSKCF_struct.maskArray,3)
                    
                    tmpMaskArray(:,:,i)=...
                        roiFromBB(shapeDSKCF_struct.maskArray(:,:,i),bbIn);
                end
                shapeDSKCF_struct.maskArray=padarray(tmpMaskArray,[segmentHIncrement segmentWIncrement]);
                
                
            end
            
            if(sizeOfSegmenter<sizeOfTarget*0.9  && minSizeOK )
                estimatedShapeBB=enlargeBB(estimatedShapeBB ,-0.05,imageSize);
                shapeDSKCF_struct.growingStatus=false;
                newOutput=true;
            end
            
            
        end
        
        
    end
end

end

