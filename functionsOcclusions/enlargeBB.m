function bb = enlargeBB(smallBB,a,size)
% ENLARGEBB.m enlarges the size of a bounding box
%
%   ENLARGEBB is a function to enlarge the size of a bounding box. Please
%   note that  this function was part of the RGBD tracker code presented in
%   [1] and available under under Open Source MIT License at
%    http://tracking.cs.princeton.edu/code.html
%
%   INPUT:
%   - smallBB   input bounding box in the format [topLeftX, topLeftY,
%   bottomRightX, bottomRightY] read as [columnIndexTopLeft, rowIndexTopLeft,
%   columnIndexBottomRight, rowIndexBottomRight]
%   -size size of the image
%   -a increasing factor (i.e 0.1 for an increase of 10 percent)
%   OUTPUT
%   - bb    output bounding box in the format [topLeftX, topLeftY,
%   bottomRightX, bottomRightY] read as [columnIndexTopLeft, rowIndexTopLeft,
%   columnIndexBottomRight, rowIndexBottomRight]
%
%  See also CHECKOCCLUSIONSDSKCF_NOISEMODEL,
%  CHECKOCCLUSIONSDSKCF_SECONDPLANE,TARGETSEARCHDSKCF, SINGLEFRAMEDSKCF
%
%  [1] Shuran Song and Jianxiong Xiao. Tracking Revisited using RGBD
%  Camera: Baseline and Benchmark. 2013.
%

x = a * (smallBB(3)-smallBB(1));
y = a * (smallBB(4)-smallBB(2));
bb(1)=max(1,smallBB(1)-x);
bb(2)=max(1,smallBB(2)-y);
bb(3)=min(size(2),smallBB(3)+x);
bb(4)=min(size(1),smallBB(4)+y);
bb=round(bb(:));
end

