function [occBB] = occludingObjectSegDSKCF(depthIm,trackerDSKCF_struct)
%OCCLUDINGOBJECTSEGDSKCF function for segmenting the occluding object
%
%OCCLUDINGOBJECTSEGDSKCF.m is the function that segments the occluding
%object. For more information about how DSKCF handles occlusions see [1].
%Please note that  this function was partially built extending the RGBD
%tracker code presented in [2] and available under under Open Source MIT
%License at
% http://tracking.cs.princeton.edu/code.html
%
%
%  INPUT:
%  - depthIm   current depth image (16BIT)
%  - trackerDSKCF_struct  DS-KCF tracker data structure
%
%
%  OUTPUT
%  - occBB Bounding box of the occluding object in the format [topLeftX,
%  topLeftY, bottomRightX, bottomRightY] read as [columnIndexTopLeft,
%   rowIndexTopLeft, columnIndexBottomRight, rowIndexBottomRight]
%
% See also BB_OVERLAP, ENLARGEBB, ROIFROMBB
%
%
%  [1] S. Hannuna, M. Camplani, J. Hall, M. Mirmehdi, D. Damen, T.
%  Burghardt, A. Paiement, L. Tao, DS-KCF: A real-time tracker for RGB-D
%  data, Journal of Real-Time Image Processing
%
%  [2] Shuran Song and Jianxiong Xiao. Tracking Revisited using RGBD
%  Camera: Baseline and Benchmark. 2013.
%
%  University of Bristol
%  Massimo Camplani and Sion Hannuna
%
%  massimo.camplani@bristol.ac.uk
%  hannuna@compsci.bristol.ac.uk

bb=enlargeBB(trackerDSKCF_struct.currentTarget.bb ,0.05,size(depthIm));
%bb=enlargeBB(bb ,0.05,size(depthMapCurr));
selectedPix=trackerDSKCF_struct.currentTarget.LabelRegions==trackerDSKCF_struct.currentTarget.regionIndex;
front_depth=roiFromBB(depthIm,bb);
tmpMean=trackerDSKCF_struct.currentTarget.Centers(trackerDSKCF_struct.currentTarget.regionIndex);%mean(depthIm(selectedPix));
depthVector=double(front_depth(selectedPix));
tmpStd=std(depthVector);
if(tmpStd<5)
    tmpStd=trackerDSKCF_struct.currentTarget.stdDepthObj;
end

occmask=abs(double(depthIm)-tmpMean)<tmpStd;
occmask= occmask & depthIm>0;
%find the main connected component
tarBBProp=regionprops(occmask,'BoundingBox','Area');
if(isempty(tarBBProp))
    occBB=[];
else if(length(tarBBProp)==1)
        occBB=tarBBProp.BoundingBox;
        bbVector=cat(1, tarBBProp.BoundingBox)';
        overlap = bb_overlap(bb,bbVector);
        %use extrema points.....
        if(overlap>0.15)
            occBB=ceil([occBB(1), occBB(2),occBB(1)+occBB(3),occBB(2)+occBB(4)]);
        else
            occBB=[];
        end
        
    else
        areas= cat(1, tarBBProp.Area);
        bbVector=cat(1, tarBBProp.BoundingBox)';
        %clean small areas....
        minArea=trackerDSKCF_struct.currentTarget.w*trackerDSKCF_struct.currentTarget.h*0.05;
        areaSmallIndex=areas<minArea;
        %exclude the small area index!!!!!!!!
        bbVector(:,areaSmallIndex) = [];
        bbVector(3:4,:)=bbVector(1:2,:)+bbVector(3:4,:);
        areas(areaSmallIndex)= [];
        overlap = bb_overlap(bb,bbVector);
        [maxV,maxIndex]=max(overlap);
        if(maxV>-100)%(maxV>0)
            occBB=bbVector(:,maxIndex)';
        else
            occBB=[];
        end
        
    end
    
    occBB=occBB';
    
    
end

