% FROMCENTRALPOINTTOBB.m is a function for calculating target bounding box
%
%  FROMCENTRALPOINTTOBB is a function that calculates the target bounding
%  box given the centroid position and the size of the target
%
%   INPUT:
%  -centerX target's centroid coordinate (column in the image plane)
%  -centerY target's centroid coordinate (row in the image plane)
%  -width target's width
%  -height target's height
% -maxX,maxY image limits
%
%   OUTPUT
%  -bb calculated bounding box in the format [topLeftY, topLeftX,
%   bottomRightY, bottomRightX] read as [rowIndexTopLeft, columnIndexTopLeft,
%   rowIndexBottomRight, columnIndexBottomRight]
%
%  See also SINGLEFRAMEDSKCF, FROMCENTRALPOINTTOBB 
%
%  University of Bristol
%  Massimo Camplani and Sion Hannuna
%  
%  massimo.camplani@bristol.ac.uk
%  hannuna@compsci.bristol.ac.uk

function bb=fromCentralPointToBB(centerX,centerY,width,height,maxX,maxY)

bb(1)=max(1,centerX-width/2);%column indexes
bb(2)=max(1,centerY-height/2);%row indexes
bb(3)=min(maxX,centerX+width/2);%column indexes
bb(4)=min(maxY,centerY+height/2);%row indexes
bb=floor(bb(:));
end