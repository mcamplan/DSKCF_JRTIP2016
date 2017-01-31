function [dsKCFoutputSr,dsKCFoutputSq,dsKCFsegmentationOut, avTime,totalTime,timeMatrix] = ...
    wrapperDSKCF(video_path, depth_path, img_files, depth_files, pos, target_sz, ...
    DSKCFparameters, show_visualization,save_Images,dest_path,noBitShift)
% WRAPPERDSKCF.m is the wrapper function for the DS-KCF tracker [1]
%
%   WRAPPERDSKCF is the wrapper function of the DS-KCF tracker [1]. This
%   function, given the tracker parameters and the information about the
%   data that will be processed, initializes the tracker's data structures
%   and starts the frame by frame tracking processing. WRAPPERDSKCF
%   eventually save/visualize the tracker image output (images with
%   overlayed tracker bounding box) in the specified folder. Please note
%   that  this function was partially built extending the KCF tracker code
%   presented by Joao F. Henriques, in http://www.isr.uc.pt/~henriques/.
%
%   INPUT:
%  -video_path absolute path where color data is stored -depth_path
%  absolute path where depth data is stored 
%  -img_files the list of Color data files (N is the number of frame to be
%  processed) -depth_files  the list of Depth data files (N is the number
%  of frame to be processed)
%  -pos initial DS-KCF tracker position pos=[y x] where x is the column
%  index and y is the row index of the image
%  -target_sz initial target size target_sz=[height,width]
%  -DSKCFparameters structure containing DSKCF parameters, see the test
%  script testDS-KCFScripts\runDSKCF.m
%  -show_visualization flag to show the tracking results live in a matlab
%  figure
%  -save_Images save tracker's output as images with overlayed tracker
%   bounding box) in the specified folder dest_path
%  -dest_path absolute path of the destination folder where tracker's
%  output is saved
%  -noBitShift boolean value to properly load the depth data. In some cases
%  the depth is stored in 16bit images where every pixel contains the depth
%  data in mm. In this case noBitShift has to be set to false. On the
%  contrary for other datasets (i.e. Princeton RGB-D) the data stored need
%  to be shifted as shown in this function
%
%   OUTPUT
%  -dsKCFoutputSr tracker's output using scale factor Sr in [1]. This a Nx5
%  vector containing in each row the tracker results for the corresponding
%  frame. In particular the output is formatted as suggested in the
%  princetonRGB-D dataset in [3]. For each row the first four columns
%  contain the bounding box information in the format [topLeftX, topLeftY,
%  bottomRightX, bottomRightY]. Note that in case of lost target the row
%  contains NaN values. The fifth column contains a flag to indicate
%  occlusions cases.
%
%  -dsKCFoutputSq tracker's output using scale factor Sq in [1]. This a Nx5
%  vector containing in each row the tracker results for the corresponding
%  frame. In particular the output is formatted as suggested in the
%  princetonRGB-D dataset in [3]. For each row the first four columns
%  contain the bounding box information in the format [topLeftX, topLeftY,
%  bottomRightX, bottomRightY]. Note that in case of lost target the row
%  contains NaN values. The fifth column contains a flag to indicate
%  occlusions cases.
%
%  -dsKCFsegmentationOut tracker's output using scale factor shape module
%  in [1]. This a Nx5 vector containing in each row the tracker results for
%  the corresponding frame. In particular the output is formatted as
%  suggested in the princetonRGB-D dataset in [3]. For each row the first
%  four columns contain the bounding box information in the format
%  [topLeftX, topLeftY, bottomRightX, bottomRightY]. Note that in case of
%  lost target the row contains NaN values. The fifth column contains a
%  flag to indicate occlusions cases.
%   
%  -avTime this is a Nx1 vector containing the processing time, expressed
%  in ms,for each frame. Note that it does not condiser the reading and
%  writing time from the disk
%
%  -totalTime is the total time required to process the entire sequence
%  
%  -timeMatrix NxM matrix containing the processing rate for all the main
%  modules of
%   the ds-kcf (see [1] for more details)
% 
%  See also LOAD_VIDEO_INFO_BOBOTRESULTS,LOAD_VIDEO_INFO_DEPTHFROMMAT,
%  INITDSKCFPARAM,INITDSKCFTRACKER,FROMCENTRALPOINTTOBB,INITDISTRIBUTIONFAST
%  MANUALBBDRAW_OCC, MANUALBBDRAW_OCC_WITHLABELSVISUALIZE,
%  INITDSKCFTRACKER_OCCLUDER,SINGLEFRAMEDSKCF
%
%
%  [1] S. Hannuna, M. Camplani, J. Hall, M. Mirmehdi, D. Damen, T.
%  Burghardt, A. Paiement, L. Tao, DS-KCF: A real-time tracker for RGB-D
%  data, Journal of Real-Time Image Processing
%
%  [2] J. F. Henriques, R. Caseiro, P. Martins, and J. Batista. High-speed
%  tracking with kernelized correlation filters. Pattern Analysis and
%  Machine Intelligence, IEEE Transactions on, 2015.
%
%  [3] Shuran Song and Jianxiong Xiao. Tracking Revisited using RGBD
%  Camera: Baseline and Benchmark. 2013.
%
%  University of Bristol
%  Massimo Camplani and Sion Hannuna
%
%  massimo.camplani@bristol.ac.uk
%  hannuna@compsci.bristol.ac.uk


%As suggested in the original code of KCF resize target and image....
resize_image = (sqrt(prod(target_sz)) >= 100);  %diagonal size >= threshold
if resize_image,
    pos = floor(pos / 2);
    target_sz = floor(target_sz / 2);
end


%window size, taking padding into account
DSKCFparameters.window_sz = floor(target_sz * (1 + DSKCFparameters.padding));

%initialize the scale parameters for the DSKCF algorithm accortind to the
%selected Sq (see [1]) information contained in DSKCFparameters
scaleDSKCF_struct=initDSKCFparam(DSKCFparameters,target_sz,pos);

%initialize shape struct
shapeDSKCF_struct=initDSKCFshape(5,0);

%check if the scale is properly initialized....
if(isempty(scaleDSKCF_struct))
    disp('Scale structure initialization failed, tracking aborted');
    dsKCFoutputSr=[];
    dsKCFoutputSq=[];
    dsKCFsegmentationOut=[];
    avTime=[];
    return;
end

%note: variables ending with 'f' are in the Fourier domain.

totalTime = 0;  %to calculate FPS

%init variables
positions = zeros(numel(img_files), 2); % where store tracker centroids
dsKCFoutputSr = zeros(numel(img_files), 2);
dsKCFoutputSq = zeros(numel(img_files), 2);
dsKCFsegmentationOut = zeros(numel(img_files), 2);

%auxiliary variables to compute the final results of the tracker, in fact
%DS-KCF tracking core returns only the centroid of the target, we need also
%to store the size of it to return proper output parameters
sizeSr = zeros(numel(img_files), 2);
sizeSq = zeros(numel(img_files), 2);

frameCurr=[]; %contains depth and color data of the current frame depth16Bit
%contains the 16bits depth information in mm. depth contains
%the normalize depth data as a grayscale image coded with 8bits
%depthNoData is the mask to identify missing depth data
%rgb is the color image and gray is the grayscale version of it

framePrev=[]; %contains the same information of frameCurr but related
%to the previous frame


occlusionState=[];% vector containing flags about the detected occlusion state
%%%% TO BE CHECKED
nanPosition=[];

avTime=[];

numberTimeIntervals=8;%7 interval plus 1 total
timeMatrix=repmat(-1,numel(img_files),numberTimeIntervals);%no es

segmentedBB=[];
segmentedSize=[];

%%FRAME BY FRAME TRACKING.....
for frame = 1:numel(img_files),
    %load images
    im = imread([video_path img_files{frame}]);
    depth = imread([depth_path depth_files{frame}]);
    
    %% inserting type control for depth image
    if(isa(depth,'uint16'))
        if(noBitShift==false)
            depth = bitor(bitshift(depth,-3), bitshift(depth,16-3));
        end
        %depth data in mm
        depth16Bit = depth;
        
        %Normalize depth data as a grayscale image [0 255]
        depth = double(depth);
        depth(depth==0) = 10000;
        depth = (depth-500)/8500;%only use the data from 0.5-8m
        depth(depth<0) = 0;
        depth(depth>1) = 1;
        depth = uint8(255*(1 - depth));
    end
    
    %resize images
    if size(im,3) > 1,
        imRGB=im;
        im = rgb2gray(im);
    else
        imRGB=im;
        imRGB(:,:,2)=im;
        imRGB(:,:,3)=im;
    end
    
    if resize_image,
        im = imresize(im, 0.5);
        imRGB = imresize(imRGB, 0.5);
        depth = imresize(depth, 0.5);
        depth16Bit = depth16Bit((1:2:end),(1:2:end));
    end
    
    
    
    %start measuring the time!!!!
    tTotal=tic();
    firstFrame=frame==1;
    
    %Insert current frame data
    frameCurr.rgb   = imRGB;
    frameCurr.gray   = im;
    frameCurr.depth = double(depth);
    frameCurr.depthNoData=depth16Bit==0;
    frameCurr.depth16Bit=depth16Bit;
    
    %for the first frame initialize the structures
    if(firstFrame)
        segmentedMASK=repmat(0,size(frameCurr.depth));
        trackerDSKCF_struct=initDSKCFtracker();
        %check if the scale is properly initialized....
        if(isempty(trackerDSKCF_struct))
            disp('DS-KCF tracker structure initialization failed, tracking aborted');
            dsKCFoutputSr=[];
            dsKCFoutputSq=[];
            avTime=[];
            return;
        end
        %%INITIALIZE HISTOGRAMS....
        framePrev.rgb   = imRGB;
        framePrev.gray   = im;
        framePrev.depth = depth;
        framePrev.depthNoData=depth16Bit==0;
        framePrev.depth16Bit=depth16Bit;
        
        trackerDSKCF_struct.previousTarget.posX=pos(2);
        trackerDSKCF_struct.previousTarget.posY=pos(1);
        trackerDSKCF_struct.previousTarget.h=scaleDSKCF_struct.target_sz(scaleDSKCF_struct.i).target_sz(1);
        trackerDSKCF_struct.previousTarget.w=scaleDSKCF_struct.target_sz(scaleDSKCF_struct.i).target_sz(2);
        trackerDSKCF_struct.previousTarget.bb=fromCentralPointToBB...
            (trackerDSKCF_struct.previousTarget.posX,trackerDSKCF_struct.previousTarget.posY,...
            trackerDSKCF_struct.previousTarget.w,trackerDSKCF_struct.previousTarget.h,size(im,2),size(im,1));
        trackerDSKCF_struct.currentTarget.meanDepthObj=0;% mean depth of the tracker object
        %initialize depth distributions
        [trackerDSKCF_struct.previousTarget.meanDepthObj,trackerDSKCF_struct.previousTarget.stdDepthObj,...
            trackerDSKCF_struct.previousTarget.LabelRegions,...
            trackerDSKCF_struct.previousTarget.regionIndex,...
            trackerDSKCF_struct.previousTarget.Centers,...
            trackerDSKCF_struct.previousTarget.LUT] = ...
            initDistributionFast(trackerDSKCF_struct.previousTarget.bb, ...
            framePrev.depth16Bit,framePrev.depthNoData);
        
        %for the first frame copy everything also in the current target
        trackerDSKCF_struct.currentTarget=trackerDSKCF_struct.previousTarget;
        
        %set the depth of the initial target in the scale data structure
        scaleDSKCF_struct.InitialDepth = trackerDSKCF_struct.previousTarget.meanDepthObj;
        scaleDSKCF_struct.currDepth = trackerDSKCF_struct.previousTarget.meanDepthObj;
        
        %initialize structures for the occluder object
        trackerDSKCF_structOccluder=initDSKCFtracker_occluder();
        DSKCFparameters_Occluder=DSKCFparameters;%these need to be resetted eventually in some parts
        
        %figure initialization
        if(show_visualization)
            
            myFigColor=figure();
            myFigDepth=figure();
            set(myFigDepth,'resize','off');
            set(myFigColor,'resize','off');
        end
        
        %take segmentation results for the first frame
        trackerDSKCF_struct.currentTarget.segmentedBB=...
            trackerDSKCF_struct.currentTarget.bb';
    end

    %DS-KCF tracker code need as input the position expressed as [y x],
    %remember this particular while reading the code!!!!!

    [pos,trackerDSKCF_struct,trackerDSKCF_structOccluder,scaleDSKCF_struct,...
        DSKCFparameters_Occluder,segmentedMASK,shapeDSKCF_struct,timeMatrix(frame,1:7)]=...
        singleFrameDSKCF(firstFrame,pos,frameCurr,trackerDSKCF_struct,DSKCFparameters,...
        scaleDSKCF_struct,trackerDSKCF_structOccluder,DSKCFparameters_Occluder,shapeDSKCF_struct);
    
    if(isempty(trackerDSKCF_struct.currentTarget.segmentedBB))
        segmentedBB=[segmentedBB;segmentedBB(end,:)];
    else
        segmentedBB=[segmentedBB;trackerDSKCF_struct.currentTarget.segmentedBB];
    end
    %Compose the Occlusion state vector for results
    occlusionState=[occlusionState ;trackerDSKCF_struct.currentTarget.underOcclusion];
    
    %tracking finished,
    avTime=[avTime; toc(tTotal)];
    timeMatrix(frame,8)=avTime(end);
    totalTime = totalTime + avTime(end);
    
    %% Just visualize......
    if (save_Images==false && show_visualization==true)
        
        %eventually re-scale the images
        if(resize_image)
            imRGB = imresize(imRGB, 2);
            depth = imresize(depth, 2);
        end
        
        %empty tracking, so mark this frame
        if(isempty(pos))
            bbToPlot=[];
            nanPosition=[nanPosition; 1];
        else
            nanPosition=[nanPosition; 0];
            
            %use the Sr scale factor (see [1] for more details)
            sr = scaleDSKCF_struct.InitialDepth / scaleDSKCF_struct.currDepth;
            
            targ_sz = round(scaleDSKCF_struct.InitialTargetSize * sr);
            
            %calculate the corresponding bounding box for Plotting!!!!
            %in this case we need [topLeftX, topLeftY,W,H]
            bbToPlot = [pos([2,1]) - targ_sz([2,1])/2, targ_sz([2,1])];
            if(resize_image)
                bbToPlot=bbToPlot*2;
            end
        end
        
        bbOCCToPlot=[];
        if(trackerDSKCF_struct.currentTarget.underOcclusion)
            widthOCC=trackerDSKCF_struct.currentTarget.occBB(3)-trackerDSKCF_struct.currentTarget.occBB(1);
            heightOCC=trackerDSKCF_struct.currentTarget.occBB(4)-trackerDSKCF_struct.currentTarget.occBB(2);
            bbOCCToPlot=[trackerDSKCF_struct.currentTarget.occBB(1:2); widthOCC; heightOCC]';
            if(resize_image)
                bbOCCToPlot=bbOCCToPlot*2;
            end
        end
        
        if(frame==1)
            manualBBdraw_OCC_WithLabelsVisualize(imRGB,bbToPlot,bbOCCToPlot,'r','y',4,'DS-KCF','Occluder',myFigColor);
            positionColor=get(gcf,'OuterPosition');
            positionColor(1)=positionColor(1)-floor(positionColor(3)/2) -25;
            set(gcf,'OuterPosition',positionColor);
            manualBBdraw_OCC_WithLabelsVisualize(depth,bbToPlot,bbOCCToPlot,'r','y',4,'DS-KCF','Occluder',myFigDepth);
            positionDepth=get(gcf,'OuterPosition');
            positionDepth(1)=positionDepth(1)+floor(positionDepth(3)/2) +25;
            set(gcf,'OuterPosition',positionDepth);
        else
            %myFigColor=figure();
            %myFigDepth=figure();
            clf(myFigColor);
            manualBBdraw_OCC_WithLabelsVisualize(imRGB,bbToPlot,bbOCCToPlot,'r','y',4,'DS-KCF','Occluder',myFigColor);
            
            clf(myFigDepth);
            manualBBdraw_OCC_WithLabelsVisualize(depth,bbToPlot,bbOCCToPlot,'r','y',4,'DS-KCF','Occluder',myFigDepth);
            drawnow
            pause(0.05)
        end
        
    end
    %% Visualize and save.....
    if (save_Images==true && show_visualization==true)
        
        %eventually re-scale the images
        if(resize_image)
            imRGB = imresize(imRGB, 2);
            depth = imresize(depth, 2);
        end
        
        %empty tracking, so mark this frame
        if(isempty(pos))
            bbToPlot=[];
            nanPosition=[nanPosition; 1];
        else
            nanPosition=[nanPosition; 0];
            
            %use the Sr scale factor (see [1] for more details)
            sr = scaleDSKCF_struct.InitialDepth / scaleDSKCF_struct.currDepth;
            
            targ_sz = round(scaleDSKCF_struct.InitialTargetSize * sr);
            
            %calculate the corresponding bounding box for Plotting!!!!
            %in this case we need [topLeftX, topLeftY,W,H]
            bbToPlot = [pos([2,1]) - targ_sz([2,1])/2, targ_sz([2,1])];
            if(resize_image)
                bbToPlot=bbToPlot*2;
            end
        end
        
        bbOCCToPlot=[];
        if(trackerDSKCF_struct.currentTarget.underOcclusion)
            widthOCC=trackerDSKCF_struct.currentTarget.occBB(3)-trackerDSKCF_struct.currentTarget.occBB(1);
            heightOCC=trackerDSKCF_struct.currentTarget.occBB(4)-trackerDSKCF_struct.currentTarget.occBB(2);
            bbOCCToPlot=[trackerDSKCF_struct.currentTarget.occBB(1:2); widthOCC; heightOCC]';
            if(resize_image)
                bbOCCToPlot=bbOCCToPlot*2;
            end
        end
        
        if(frame==1)
            imRGB_tracked=manualBBdraw_OCC_WithLabelsVisualize(imRGB,bbToPlot,bbOCCToPlot,'r','y',4,'DS-KCF','Occluder',myFigColor);
            positionColor=get(gcf,'OuterPosition');
            positionColor(1)=positionColor(1)-floor(positionColor(3)/2) -25;
            set(gcf,'OuterPosition',positionColor);
            dept_tracked=manualBBdraw_OCC_WithLabelsVisualize(depth,bbToPlot,bbOCCToPlot,'r','y',4,'DS-KCF','Occluder',myFigDepth);
            positionDepth=get(gcf,'OuterPosition');
            positionDepth(1)=positionDepth(1)+floor(positionDepth(3)/2) +25;
            set(gcf,'OuterPosition',positionDepth);
        else
            %myFigColor=figure();
            %myFigDepth=figure();
            clf(myFigColor);
            imRGB_tracked=manualBBdraw_OCC_WithLabelsVisualize(imRGB,bbToPlot,bbOCCToPlot,'r','y',4,'DS-KCF','Occluder',myFigColor);
            
            clf(myFigDepth);
            dept_tracked=manualBBdraw_OCC_WithLabelsVisualize(depth,bbToPlot,bbOCCToPlot,'r','y',4,'DS-KCF','Occluder',myFigDepth);
            drawnow
            pause(0.05)
        end
        
        tmpColorName=[dest_path '/trackedColor_' num2str(frame) '.jpg'];
        if(frame==1)
            imRGB_trackedMASK=imRGB_tracked(:,:,1)==204;
            imRGB_trackedMASK=imRGB_trackedMASK & imRGB_tracked(:,:,2)==204;
            imRGB_trackedMASK=imRGB_trackedMASK & imRGB_tracked(:,:,3)==204;
            finalMask=imfill(~imRGB_trackedMASK,'holes');
            rp=regionprops(finalMask,'boundingbox','Area');
            
            areaList= cat(1, rp.Area);
            bbList= cat(1, rp.BoundingBox);
            [ddd,areaMax]=max(areaList);
            rect=floor(bbList(areaMax,:));
            
        end
        %clean the image from gray fig border....and save
        imRGB_trackedNEW=imRGB_tracked(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3),:);
        imwrite(imRGB_trackedNEW,tmpColorName);
        
        %clean the image from gray fig border.... and save
        tmpDepthName=[dest_path '/trackedDepth_' num2str(frame) '.jpg'];
        dept_trackedNEW=dept_tracked(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3),:);
        imwrite(dept_trackedNEW,tmpDepthName);
    end
    %% just save images
    if (save_Images==true && show_visualization==false)
        
        %eventually re-scale the images
        if(resize_image)
            imRGB = imresize(imRGB, 2);
            depth = imresize(depth, 2);
            segmentedMASK= imresize(segmentedMASK,2);
        end
        
        tmpMaskName=[dest_path '/binaryMask_' num2str(frame) '.jpg'];
        imwrite(segmentedMASK>0,tmpMaskName);
        
        tmpMaskName=[dest_path '/binaryMaskColored_' num2str(frame) '.jpg'];
        segmentedMASKColored = imRGB.*(uint8(repmat(segmentedMASK>0,[1,1,3])));
        imwrite(segmentedMASKColored,tmpMaskName);
        
        bbSegToPlot=[1 1 1 1];
        if(isempty(trackerDSKCF_struct.currentTarget.segmentedBB)==false)
            widthOCC=trackerDSKCF_struct.currentTarget.segmentedBB(3)-trackerDSKCF_struct.currentTarget.segmentedBB(1);
            heightOCC=trackerDSKCF_struct.currentTarget.segmentedBB(4)-trackerDSKCF_struct.currentTarget.segmentedBB(2);
            bbSegToPlot=[trackerDSKCF_struct.currentTarget.segmentedBB(1:2)'; widthOCC; heightOCC]';
            if(resize_image)
                bbSegToPlot=bbSegToPlot*2;
            end
        end
        
        
        %empty tracking, so mark this frame
        if(isempty(pos))
            bbToPlot=[1 1 1 1];
            nanPosition=[nanPosition; 1];
        else
            nanPosition=[nanPosition; 0];
            
            %use the Sr scale factor (see [1] for more details)
            sr = scaleDSKCF_struct.InitialDepth / scaleDSKCF_struct.currDepth;
            
            targ_sz = round(scaleDSKCF_struct.InitialTargetSize * sr);
            
            %calculate the corresponding bounding box for Plotting!!!!
            %in this case we need [topLeftX, topLeftY,W,H]
            bbToPlot = [pos([2,1]) - targ_sz([2,1])/2, targ_sz([2,1])];
            %bbToPlot=[segmentedBB(end,[1 2]),(segmentedBB(end,[3 4])-segmentedBB(end,[1 2]))];
            if(resize_image)
                bbToPlot=bbToPlot*2;
            end
        end
        
        bbOCCToPlot=[1 1 1 1];
        if(trackerDSKCF_struct.currentTarget.underOcclusion)
            widthOCC=trackerDSKCF_struct.currentTarget.occBB(3)-trackerDSKCF_struct.currentTarget.occBB(1);
            heightOCC=trackerDSKCF_struct.currentTarget.occBB(4)-trackerDSKCF_struct.currentTarget.occBB(2);
            bbOCCToPlot=[trackerDSKCF_struct.currentTarget.occBB(1:2); widthOCC; heightOCC]';
            if(resize_image)
                bbOCCToPlot=bbOCCToPlot*2;
            end
        end
        
        bbSeparate(1,:)=bbToPlot;
        bbSeparate(2,:)=bbOCCToPlot;
        bbSeparate(3,:)=bbSegToPlot;
        
        colorVector='ryb';
        myText{1}='DS-KCF';
        myText{2}='Occluder';
        myText{3}='DS-KCF+SEG';
        myText={};
        myTextColor='ryb';
        
        %draw the bb square....
        imRGB_tracked=manualBBdraw_MULTIPLE(imRGB,bbSeparate,colorVector,4,myText,myTextColor);
        dept_tracked=manualBBdraw_MULTIPLE(depth,bbSeparate,colorVector,4,myText,myTextColor);

        %end
        tmpColorName=[dest_path '/trackedColor_' num2str(frame) '.jpg'];
        if(frame==1)
            imRGB_trackedMASK=imRGB_tracked(:,:,1)==204;
            imRGB_trackedMASK=imRGB_trackedMASK & imRGB_tracked(:,:,2)==204;
            imRGB_trackedMASK=imRGB_trackedMASK & imRGB_tracked(:,:,3)==204;
            finalMask=imfill(~imRGB_trackedMASK,'holes');
            rp=regionprops(finalMask,'boundingbox','Area');
            
            areaList= cat(1, rp.Area);
            bbList= cat(1, rp.BoundingBox);
            [ddd,areaMax]=max(areaList);
            rect=floor(bbList(areaMax,:));
            
        end
        %clean the image from gray fig border....and save
        imRGB_trackedNEW=imRGB_tracked(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3),:);
        imwrite(imRGB_trackedNEW,tmpColorName);
        
        %clean the image from gray fig border.... and save
        tmpDepthName=[dest_path '/trackedDepth_' num2str(frame) '.jpg'];
        dept_trackedNEW=dept_tracked(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3),:);
        imwrite(dept_trackedNEW,tmpDepthName);
    end
    
    %now generate the results, starting from the tracker output!!!
    % the object has being tracked....
    if(isempty(pos)==false)
        %accumulate the position of the DS-KCF tracker remember format [y x]
        positions(frame,:) = pos;
        %take the sq size (see [1]) invert the coordinate as you need to
        %combine this with positions
        sqVector(frame,:)= scaleDSKCF_struct.target_sz(scaleDSKCF_struct.i).target_sz([2,1]);
        
        %use the Sr scale factor (see [1] for more details)
        sr = scaleDSKCF_struct.InitialDepth / scaleDSKCF_struct.currDepth;
        targ_sz = round(scaleDSKCF_struct.InitialTargetSize * sr);
        %invert the coordinate as you need to combine this with positions
        srVector(frame,:) = targ_sz([2,1]);
        
    else
        pos=positions(frame-1,:);
        positions(frame,:) = pos;
        
        tmpSize=sqVector(frame-1,:);
        sqVector(frame,:)= tmpSize;
        
        tmpSize=srVector(frame-1,:);
        srVector(frame,:) = tmpSize;
        
    end
    
    %Update PAST Target data structure
    if(frame>1)
        %previous target entries
        trackerDSKCF_struct.previousTarget=trackerDSKCF_struct.currentTarget;
    end
    
    
end

if resize_image,
    positions = positions * 2;
    sqVector= sqVector*2;
    srVector= srVector*2;
    segmentedBB=2*segmentedBB;
end

%now generate the final results, this are in the format requested by the
%Princeton RGB-D dataset in the format [topLeftX, topLeftY, bottomRightX,
%bottomRightY]. Then move from the change from the center + target size
%format to the above mentioned one

dsKCFoutputSr=[positions(:,[2,1]) - srVector/2, positions(:,[2,1]) + srVector/2];
dsKCFoutputSq=[positions(:,[2,1]) - sqVector/2, positions(:,[2,1]) + sqVector/2];
dsKCFsegmentationOut=segmentedBB;

%now set to NaN the output for the frames where the tracker was not
%available
dsKCFoutputSr(nanPosition>0,:)=NaN;
dsKCFoutputSq(nanPosition>0,:)=NaN;
dsKCFsegmentationOut(nanPosition>0,:)=NaN;

%add the occlusion state vector.
dsKCFoutputSr=[dsKCFoutputSr, occlusionState];
dsKCFoutputSq=[dsKCFoutputSq, occlusionState];
dsKCFsegmentationOut=[dsKCFsegmentationOut, occlusionState];
end

