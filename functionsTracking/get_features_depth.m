function [x separateDepth] = get_features_depth(im, depth, features, cell_size, cos_window)
%GET_FEATURES_DEPTH Extracts dense features from image
%
%GET_FEATURES_DEPTH.m is a function used for extracting dense features from
%depth and color image according to the feature used. For more information
%about the DS-KCF and KCF model update see [1,2]. Please note that  this
%function was partially built extending the KCF tracker code presented by
%Joao F. Henriques, in http://www.isr.uc.pt/~henriques/.
%
%
%  INPUT:
%  -im color data
%  -depth depth data
%  -features struct containing feature info.
%  -cell_size HOG parameter
%  -cos_window cosine window to smooth data in the Fourier domain
%
%  OUTPUT:
%  -x features extracted 
%  -separateDepth features extracted for depth in case color and depth
%  features are selected, in the other cases this matrix is empty
%
%  See also FHOG, MAXRESPONSEDSKCF, MAXRESPONSEDEPTHWEIGHTDSKCF,
%  MODELUPDATEDSKCF
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

separateDepth=[];

if features.rawDepth,
    x = double(im) ;%x = double(im) / 255;
    
    x = (x - mean(x(:)))/std(x(:));
    
end

if features.rawColor,
    
    x = double(im)/255 ;%x = double(im) / 255;
    x = x - mean(x(:));
    
end

if features.rawConcatenate,
    
    x = double(im) ;

    x = (x - mean(x(:)))/std(x(:));
    
    xD = double(depth);
    
    xD = (xD - mean(xD(:)))/std(xD(:));
    
    x = cat(3, x, xD);
end

if features.rawLinear,

    x = double(im) ;

   x = (x - mean(x(:)))/std(x(:));
   
   separateDepth = double(depth);%separateDepth = double(depth) / 255;
   %separateDepth = separateDepth - mean(separateDepth(:));
   separateDepth = (separateDepth - mean(separateDepth(:)))/std(separateDepth(:));
   %x = cat(3, x, xD);
end


if features.hog_color,
    %HOG features, from Piotr's Toolbox
    x = double(fhog(single(im) / 255, cell_size, features.hog_orientations));
    x(:,:,end) = [];  %remove all-zeros channel ("truncation feature")
   % fprintf('In Hog\n');
end

if features.hog_depth,
    %HOG features, from Piotr's Toolbox
    x = double(fhog(single(depth) / 255, cell_size, features.hog_orientations));
    x(:,:,end) = [];  %remove all-zeros channel ("truncation feature")
   % fprintf('In Hog\n');
end

if features.hog_concatenate,
    %HOG features, from Piotr's Toolbox
    x = double(fhog(single(im) / 255, cell_size, features.hog_orientations));
    x(:,:,end) = [];  %remove all-zeros channel ("truncation feature")
   % fprintf('In Hog\n');
    xD = double(fhog(single(depth) / 255, cell_size, features.hog_orientations));
    xD(:,:,end) = [];  %remove all-zeros channel ("truncation feature")
    x = cat(3, x, xD);
end

if features.hog_linear,
    %HOG features, from Piotr's Toolbox
    x = double(fhog(single(im) / 255, cell_size, features.hog_orientations));
    x(:,:,end) = [];  %remove all-zeros channel ("truncation feature")
   % fprintf('In Hog\n');
    separateDepth = double(fhog(single(depth) / 255, cell_size, features.hog_orientations));
    separateDepth(:,:,end) = [];  %remove all-zeros channel ("truncation feature")
end


%process with cosine window if needed
if ~isempty(cos_window),
    x = bsxfun(@times, x, cos_window);
end



end
