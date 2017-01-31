function [tarBB, occBB, tarlist, id,occmask] = targetSearchDSKCF(bb,...
    trackerDSKCF_struct, DSKCFparameters,im,depth,depth16Bit,...
    scaleDSKCF_struct,confValue)
%TARGETSEARCHDSKCF function for segmenting the occluding object
%
%TARGETSEARCHDSKCF.m this function segments the occluding area and find
%meaningful target candidates. For more information about how DSKCF handles
%occlusions see [1]. Please note that  this function was partially built
%extending the RGBD tracker code presented in [2] and available under under
%Open Source MIT License at
% http://tracking.cs.princeton.edu/code.html
%
%
%  INPUT:
%  - depthIm   current depth image (16BIT)
%  - trackerDSKCF_struct  DS-KCF tracker data structure (see INITDSKCFTRACKER)
%
%
%  OUTPUT
%  - occBB Bounding box of the occluding object in the format [topLeftX,
%  topLeftY, bottomRightX, bottomRightY] read as [columnIndexTopLeft,
%   rowIndexTopLeft, columnIndexBottomRight, rowIndexBottomRight]
%
% See also ENLARGEBB, REGIONSEGMENTSFAST, GET_SUBWINDOW, MAXRESPONSEDSKCF,
% SINGLEFRAMEDSKCF
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
tarBB = [];
occBB=[];
id=[];

bbIn=bb;
bb=enlargeBB(bb ,0.05,size(im));

%%take Bounding box and mask from the previous segmentation!!!!!
occmask=trackerDSKCF_struct.currentTarget.LabelRegions==trackerDSKCF_struct.currentTarget.regionIndex;
centersEstimated=trackerDSKCF_struct.currentTarget.Centers;
tarBBProp=regionprops(occmask,'BoundingBox','Area');

%one single big region....just take the bounding box
if(length(tarBBProp)==1)
    occBB=tarBBProp.BoundingBox;
    %use extrema points.....
    occBB=ceil([occBB(1), occBB(2),occBB(1)+occBB(3),occBB(2)+occBB(4)]);
    %and recenter to the entire image coordinate
    occBB([1 3])=occBB([1 3])+bb(1);
    occBB([2 4])=occBB([2 4])+bb(2);
    %else select the biggest one
elseif(length(tarBBProp)>1)
    areas= cat(1, tarBBProp.Area);
    [maxV,maxIndex]=max(areas);
    occBB=tarBBProp(maxIndex).BoundingBox;
    occBB=ceil([occBB(1), occBB(2),occBB(1)+occBB(3),occBB(2)+occBB(4)]);
    
    occBB([1 3])=occBB([1 3])+bb(1);
    occBB([2 4])=occBB([2 4])+bb(2);
    
else
    %EMPTY AREA....
    occBB=[];
end

occBB=occBB';

%retrieve the list of possible candidates according to the segmentation
[tarBBList, areaList] = regionSegmentsFast( trackerDSKCF_struct.currentTarget.LabelRegions,bb);

%exclude the occluding index!!!!!!!!
if(trackerDSKCF_struct.currentTarget.regionIndex~=6666666)
    tarBBList(:,trackerDSKCF_struct.currentTarget.regionIndex) = [];
    areaList(trackerDSKCF_struct.currentTarget.regionIndex)= [];
    centersEstimated(trackerDSKCF_struct.currentTarget.regionIndex) =[];
end

minArea=trackerDSKCF_struct.currentTarget.w*...
    trackerDSKCF_struct.currentTarget.h*0.05;
areaSmallIndex=areaList<minArea;
%exclude the small area index!!!!!!!!
tarBBList(:,areaSmallIndex) = [];
areaList(areaSmallIndex)= [];
centersEstimated(areaSmallIndex)=[];

%for each target bb caculate confidience
tarlist = struct('bb',[],'Conf_color',[],'Conf_class',[]);
num_tar = size(tarBBList,2);
if isempty(tarBBList),
    %disp('total occ');
else
    tarlist.bb=nan(4,num_tar);
    tarlist.Conf_class=nan(1,num_tar);
    tarlist.Area=nan(1,num_tar);
    
    smoothFactorD_vector=[];
    for j=1:size(tarBBList,2),
        %conf  = hogGetConf(tarBBList(:,j), featurePym, svm);
        tmpBB=tarBBList(:,j);
        tmpWidthTarget=tmpBB(3)-tmpBB(1)+1;
        tmpHeightTarget=tmpBB(4)-tmpBB(2)+1;
        tmpCenter=ceil([tmpBB(1)+tmpWidthTarget/2 tmpBB(2)+tmpHeightTarget/2]);
        
        %take image patch....
        patch = get_subwindow(im, tmpCenter(2:-1:1), scaleDSKCF_struct.windows_sizes(scaleDSKCF_struct.i).window_sz);
        patch_depth = get_subwindow(depth, tmpCenter(2:-1:1), scaleDSKCF_struct.windows_sizes(scaleDSKCF_struct.i).window_sz);
        
        [ response, maxResponse, maxPositionImagePlane] = maxResponseDSKCF...
            ( patch,patch_depth, DSKCFparameters.features,DSKCFparameters.kernel,...
            tmpCenter(2:-1:1),DSKCFparameters.cell_size, ...
            scaleDSKCF_struct.cos_windows(scaleDSKCF_struct.i).cos_window,...
            trackerDSKCF_struct.model_xf,trackerDSKCF_struct.model_alphaf,...
            trackerDSKCF_struct.model_xDf,trackerDSKCF_struct.model_alphaDf,...
            size(im,1),size(im,2));
        
        %now maxPositionImagePlane has row index and column index so.....
        smoothFactorD=weightDistanceLogisticOnDepth(trackerDSKCF_struct.currentTarget.meanDepthObj,...
            double(depth16Bit(maxPositionImagePlane(1),maxPositionImagePlane(2))),...
            trackerDSKCF_struct.currentTarget.stdDepthObj);
        smoothFactorD_vector=[smoothFactorD_vector,smoothFactorD];
        
        conf=maxResponse*smoothFactorD;
        
        reEstimatedBB(1:2)=maxPositionImagePlane(2:-1:1) - [tmpWidthTarget,tmpHeightTarget]/2;
        reEstimatedBB(3:4)=maxPositionImagePlane(2:-1:1) + [tmpWidthTarget,tmpHeightTarget]/2;
        
        tarlist.bb(:,j)=reEstimatedBB;
        tarlist.Conf_class(j)= conf;
        
    end
    %output the most possible target BB
    idx =  tarlist.Conf_class>-999;
    tarlist.bb=tarlist.bb(:,idx);
    tarlist.Conf_class=tarlist.Conf_class(idx);
    
    if sum(idx(:))>0,
        [conf,id]=max(tarlist.Conf_class);
        if(smoothFactorD_vector(id)~=0)
            conf=conf/smoothFactorD_vector(id);
        end
        tarlist.Conf_class=tarlist.Conf_class./smoothFactorD_vector;
        tarlist.Conf_class(isnan(tarlist.Conf_class))=0;
        %if ~isempty(id)&&conf>0.3*svm.thr, tarBB=tarlist.bb(:,id);end
        if ~isempty(id)&&conf>confValue, tarBB=tarlist.bb(:,id);end
    end
end

end

function res=sigmFunction(x,A,K,Q,ni,B,M)

res=A+(K-A)./((1+Q*exp(-B*(x-M))).^(1/ni));

end

function smoothFactor=weightDistanceLogisticOnDepth(targetDepth,candidateDepth,targetSTD)

%smoothFactor=1;
dist=abs((targetDepth-candidateDepth))/(3*targetSTD);

Q=1;
ni=0.5;
B=3.2;
M=1.94;
%              sigmFunction(x,A,K,Q,ni,B,M)
smoothFactor=(1-sigmFunction(dist,0,1,Q,ni,B,M));

end