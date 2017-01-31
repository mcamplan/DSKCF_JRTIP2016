function [img_files, depth_files, pos, target_sz, ground_truth, video_path, depth_path] = load_video_info_BobotResults(base_path, video)
% LOAD_VIDEO_INFO_BOBOTRESULTS.m is automatically generates the relevant information for the video in the given path
%
%   LOAD_VIDEO_INFO_BOBOTRESULTS loads all the relevant information for the
%   sequence depth and video data in the given path for the BOBOT-D RGB-D
%   dataset [1]. This function have been created extending
%   load_video_info of the KCF matlab library presented by Joao F.
%   Henriques, in http://www.isr.uc.pt/~henriques/. A similar function need
%   to be implemented if another dataset that uses another naming
%   convention is processed. See for example LOAD_VIDEO_INFO_DEPTHFROMMAT
%   for the Princeton RGB-D dataset [2]
%
%   INPUT:
%  -base_path name of the top folder where sequences are stored
%  -video name of the sequence that will be processed
%
%   OUTPUT
%  -img_files the list of Color data files (images in the folder are
%  supposed to be *.png)
%  -depth_files the list of Depth data files (images in the folder are
%  supposed to be  16 bit *.png)
%  -pos initial DS-KCF tracker position pos=[y x] where x is the column
%  index and y is the row index of the image
%  -target_sz initial target size target_sz=[height,width]
%  -ground_truth ground truth information
%  -video_path  absolute path of color data
%  -depth_path  absolute path of depth data
%
%  See also LOAD_VIDEO_INFO_DEPTHFROMMAT
%
%  [1]Germán Martín García, Dominik A. Klein, Jörg Stückler, Simone
%  Frintrop, and Armin B. Cremers DAGM/OAGM Conference, August 28-31, 2012,
%  Graz, Austria
%
%  [2] S. Song and J. Xiao. Tracking revisited using RGBD camera: Unified
%  benchmark and baselines. In Computer Vision (ICCV), 2013 IEEE
%  International Conference on, pages 233–240, 2013.
%
%
%  University of Bristol
%  Massimo Camplani and Sion Hannuna
%
%  massimo.camplani@bristol.ac.uk
%  hannuna@compsci.bristol.ac.uk

%full path to the video's files
if base_path(end) ~= '/' && base_path(end) ~= '\',
    base_path(end+1) = '/';
end
video_path = [base_path video '/'];

filename = [video_path 'init.txt'];
f = fopen(filename);
assert(f ~= -1, ['No initial position or ground truth to load ("' filename '").'])

%the format is [x, y, width, height]
try
    ground_truth = textscan(f, '%f,%f,%f,%f', 'ReturnOnError',false);
catch  %#ok, try different format (no commas)
    frewind(f);
    ground_truth = textscan(f, '%f %f %f %f');
end
ground_truth = cat(2, ground_truth{:});
fclose(f);


ground_truth=floor(ground_truth);
%set initial position and size
target_sz = [ground_truth(1,4), ground_truth(1,3)];
pos = [ground_truth(1,2), ground_truth(1,1)] + floor(target_sz/2);

if size(ground_truth,1) == 1,
    %we have ground truth for the first frame only (initial position)
    ground_truth = [];
else
    %store positions instead of boxes
    ground_truth = ground_truth(:,[2,1]) + ground_truth(:,[4,3]) / 2;
end

filenameGTlist = [video_path video '.txt'];
tmpGT=load(filenameGTlist);
numOfFrames = size(tmpGT,1);

%from now on, work in the subfolder where all the images are


depth_path = [video_path 'depth/'];
video_path = [video_path 'rgb/'];

%general case, just list all images
img_files_ = dir([video_path '*.png']);
assert(~isempty(img_files_), 'No image files to load.')

depth_files_ = dir([depth_path '*.png']);
assert(~isempty(depth_files_), 'No depth files to load.')





for i = 0:numOfFrames-1
    img_files{i+1} = sprintf('%d_color.png',i);
    depth_files{i+1} = sprintf('%d_depth.png',i);
    
end


end

