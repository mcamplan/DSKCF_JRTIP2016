%REGIONSEGMENTSFAST function for extracting target candidates from the occluded area
%
%REGIONSEGMENTSFAST.m this function analyzes the segmented occluding area
%and extract meaningful target candidates. For more information about how
%DSKCF handles occlusions see [1].
%
%
%  INPUT:
%  - labelMatrix   image containing pixels' label corresponding to the
%  segmentation
%  - imageCoordinateOffset  bounding box containing the occluding area n
%  the format [topLeftX, topLeftY, bottomRightX, bottomRightY] read as
%  [columnIndexTopLeft, rowIndexTopLeft, columnIndexBottomRight,
%  rowIndexBottomRight]
%
%  OUTPUT - tarListFull Matrix that contains (in each column) target
%  candidates' bounding box in the format [topLeftX, topLeftY,
%  bottomRightX, bottomRightY] read as [columnIndexTopLeft,
%   rowIndexTopLeft, columnIndexBottomRight, rowIndexBottomRight]
%
% See also TARGETSEARCHDSKCF
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

function [ tarListFull, areaVector ] = regionSegmentsFast(labelMatrix, imageCoordinateOffset)
tarListFull=[];


tarBBProp=regionprops(labelMatrix,'BoundingBox','Area');
areaVector=cat(1, tarBBProp.Area);

for i=1:length(areaVector)
    
    tmpBB=tarBBProp(i).BoundingBox;
    %use extrema points.....
    tmpBB=ceil([tmpBB(1), tmpBB(2),tmpBB(1)+tmpBB(3),tmpBB(2)+tmpBB(4)]);
    %and recenter to the entire image coordinate
    tmpBB([1 3])=tmpBB([1 3])+imageCoordinateOffset(1);
    tmpBB([2 4])=tmpBB([2 4])+imageCoordinateOffset(2);
    
    tarListFull=[tarListFull, tmpBB'];
end

end


