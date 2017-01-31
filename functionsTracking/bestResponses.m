function [maxResponse, maxPositionImagePlane]=bestResponses(depth16Bit,...
    response,poolSize,cell_size,previousPos,meanDepthObj,stdDepthObj)
% BESTRESPONSES.m select the maximum response of DSKCF tracker 
% 
% 
%   BESTRESPONSES is is a function used for selecting the maximum DSKCF
%   response. The response of a pool of candidate positions is weighted
%   considering the depth information of the target. For more information
%   about the DSKCF response see [1]. 
%   Please note that  this function was partially built extending the KCF
%   tracker code presented by Joao F. Henriques, in
%   http://www.isr.uc.pt/~henriques/.
%
%   INPUT:
%  - depth16Bit   depth image
%  - response of the DSKCF tracker
%  -cell_size HOG parameter
%  -meanDepthObj mean depth value of the target object
%  -stdDepthObj depth standard deviation of the target object
%  -previousPos is the position of the tracked target in the previous frame
%  It is in the format [y, x] (read also as [rowIndex, columIndex])
%
%   OUTPUT
%   -maxResponse maximum value of the DSKCF response
%   -maxPositionImagePlane vector containing the position in the image
%   plane of the target's centroid. It is in the format [y, x] (read also
%   as [rowIndex, columIndex])
%  
%  See also MAXRESPONSEDEPTHWEIGHTDSKCF, wrapperDSKCF, initDSKCFparam,
%  INITDSKCFPARAM, INITDSKCFTRACKER, MODELUPDATEDSKCF
%
%[1] S. Hannuna, M. Camplani, J. Hall, M. Mirmehdi, D. Damen, T. Burghardt,
% A.Paiement, L. Tao, DS-KCF: A ~real-time tracker for RGB-D data, Journal
% of Real-Time Image Processing
%
%
%  University of Bristol 
%  Massimo Camplani and Sion Hannuna
%  
%  massimo.camplani@bristol.ac.uk 
%  hannuna@compsci.bristol.ac.uk

vert_deltaVect=[];
horiz_deltaVect=[];
maxPositionImagePlane=[];
%imagePlaneYVect=[];
responseVect=[];
depthVect=[];
smoothFactorD=[];

nRows=size(depth16Bit,1);
nCols=size(depth16Bit,2);

%analyze the firs #poolsize candidates
for i=1:poolSize

    [responseVect(i),tmpIndex]=max(response(:));
    
    [vert_deltaVect(i), horiz_deltaVect(i)] = find(response == responseVect(i), 1);
    %%clean for the next response....
    response(tmpIndex)=nan;
    
    %maxPosition=[vert_delta, horiz_delta];
    
    if vert_deltaVect(i) > size(response,1) / 2,  %wrap around to negative half-space of vertical axis
        vert_deltaVect(i) = vert_deltaVect(i) - size(response,1);
    end
    if horiz_deltaVect(i) > size(response,2) / 2,  %same for horizontal axis
        horiz_deltaVect(i) = horiz_deltaVect(i) - size(response,2);
    end

    maxPositionImagePlaneVector(i,:) = previousPos + cell_size * [vert_deltaVect(i) - 1, horiz_deltaVect(i) - 1];
    maxPositionImagePlaneVector(i,maxPositionImagePlaneVector(i,:)<1)=1;
    %maxPositionImagePlaneVector(i,maxPositionImagePlaneVector(i,1)>nRows)=nRows;
    if(maxPositionImagePlaneVector(i,1)>nRows)
        maxPositionImagePlaneVector(i,1)=nRows;
    end
    %maxPositionImagePlaneVector(i,maxPositionImagePlaneVector(i,2)>nCols)=nCols;
    if(maxPositionImagePlaneVector(i,2)>nCols)
        maxPositionImagePlaneVector(i,2)=nCols;
    end
    depthVect(i)=depth16Bit(maxPositionImagePlaneVector(i,1),maxPositionImagePlaneVector(i,2));
    
    % Now you have the position take the depth vector and calculate the weight
    smoothFactorD(i)=weightDistanceLogisticOnDepth(meanDepthObj,double(depthVect(i)),stdDepthObj);
end

responseVect=responseVect.*smoothFactorD;
[maxResponse, maxIndex]=max(responseVect);
maxPositionImagePlane=maxPositionImagePlaneVector(maxIndex,:);

end

%%FUNCTION FOR ADD THE DEPTH WEIGHT...
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

%SIGMOIDAL FUNCTION
function res=sigmFunction(x,A,K,Q,ni,B,M)

res=A+(K-A)./((1+Q*exp(-B*(x-M))).^(1/ni));

end