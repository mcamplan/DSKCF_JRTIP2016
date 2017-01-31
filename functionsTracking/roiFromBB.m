% ROIFROMBB.m is a function for extracting roi from images given a bounding box
%
%  ROIFROMBB is a function that extracts image roi from an image given a
%  bounding box
%
%   INPUT:
%  -imIN input image (can be grayscale, depth or color image)
%  -bbIn input bounding box in the format [topLeftY, topLeftX,
%   bottomRightY, bottomRightX] read as [rowIndexTopLeft, columnIndexTopLeft,
%   rowIndexBottomRight, columnIndexBottomRight]
%
%   OUTPUT
%  -imOUT selected ROI
%
% See also INITDISTRIBUTIONFAST, CHECKOCCLUSIONSDSKCF_NOISEMODEL,
% CHECKOCCLUSIONSDSKCF_SECONDPLANE, OCCLUDINGOBJECTSEGDSKCF,
% SINGLEFRAMEDSKCF,
%
%  University of Bristol
%  Massimo Camplani and Sion Hannuna
%  
%  massimo.camplani@bristol.ac.uk
%  hannuna@compsci.bristol.ac.uk

function imgOUT=roiFromBB(imgIN,bbIn)

imgOUT=imgIN(bbIn(2):bbIn(4),bbIn(1):bbIn(3),:);
end