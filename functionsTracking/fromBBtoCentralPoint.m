% FROMBBTOCENTRALPOINT.m is a function for calculating target centroid and size
%
%  FROMBBTOCENTRALPOINT is a function for calculating target centroid and
%  size given the bounding box of the target
%
%   OUTPUT:
%  -centerX target's centroid coordinate (column in the image plane)
%  -centerY target's centroid coordinate (row in the image plane)
%  -width target's width
%  -height target's height
%
%   INPUT
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

function [centerX,centerY,width,height]=fromBBtoCentralPoint(bb)

width=bb(3)-bb(1);
height=bb(4)-bb(2);
centerX=floor(bb(1)+width/2);%column indexes
centerY=floor(bb(2)+height/2);%row indexes

end