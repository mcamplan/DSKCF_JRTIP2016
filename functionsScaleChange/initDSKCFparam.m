function scaleDSKCF_struct=initDSKCFparam(DSKCFparameters,target_sz,pos)
% INITDSKCFPARAM.m initializes the scale data structure for DS-KCF tracker [1]
% 
%   INITDSKCFPARAM function initializes the scale data structure of the
%   DS-KCF tracker. In particular, some matrices are precomputed at the
%   different scales and other flags are intialized
%
%   INPUT:
%  -DSKCFparameters DS-KCF algorithm parameters 
%  -target_sz initial target size of the tracked object
%  -pos initial target position of the tracked object
%   OUTPUT
%  -scaleDSKCF_struct data structure that contains scales parameters and
%  precomputed matrices. In particular the field of the struct are 
%
%  + i current scale among Sq in [1]
%  + iPrev previous scale among Sq in [1]
%  + minStep minimum interval between the scales in Sq in [1]
%  + scales is the Sq vector in [1]
%  + step minimum interval between the scales in Sq (the same as minStep)
%  + updated flag set to 1 (or 0) when a change of scale is (or not)
%  required
%  + currDepth is the ratio between the initial depth target and the
%  current one. Is Sr in [1]
%  + InitialDepth initial ratio
%  + InitialTargetSize size of the target in the initial frame
%  + windows_sizes vector containing precomputed windows size (according
%  the padding parameter) for each scale
%  + target_sz vector containing the expected target size at the different
%  scales
%  + pos vector where is stored per each scale the target centroid position
%  + output_sigmas vector where sigma parameter is stored for each scale
%  + yfs regression targets for all the scales
%  + cos_windows precomputed cosine windows for each scale
%  + len contains the area of the target at each scale
%  + ind 
%
%  [1] S. Hannuna, M. Camplani, J. Hall, M. Mirmehdi, D. Damen, T.
%  Burghardt, A.Paiement, L. Tao, DS-KCF: A real-time tracker for RGB-D
%  data, Journal of Real-Time Image Processing
%
%
%  University of Bristol 
%  Massimo Camplani and Sion Hannuna
%  
%  massimo.camplani@bristol.ac.uk 
%  hannuna@compsci.bristol.ac.uk

scaleDSKCF_struct=[];

% Find initial scale
scaleDSKCF_struct.i = find(DSKCFparameters.scales == 1);
scaleDSKCF_struct.iPrev = scaleDSKCF_struct.i; % for inerpolating model
scaleDSKCF_struct.minStep = min(abs(diff(DSKCFparameters.scales)));
scaleDSKCF_struct.scales = DSKCFparameters.scales;
scaleDSKCF_struct.step = min(diff(DSKCFparameters.scales)); % use smallest step to decide whether to look at other scales
scaleDSKCF_struct.updated = 0;
scaleDSKCF_struct.currDepth = 1;
scaleDSKCF_struct.InitialDepth = 1;
scaleDSKCF_struct.InitialTargetSize = target_sz;

for i=1:length(DSKCFparameters.scales)
    scaleDSKCF_struct.windows_sizes(i).window_sz = round(DSKCFparameters.window_sz * DSKCFparameters.scales(i));
    scaleDSKCF_struct.target_sz(i).target_sz = round(target_sz * DSKCFparameters.scales(i));
    scaleDSKCF_struct.pos(i).pos = pos;
    
    %create regression labels, gaussian shaped, with a bandwidth
    %proportional to target size
    scaleDSKCF_struct.output_sigmas(i).output_sigma = sqrt(prod(scaleDSKCF_struct.target_sz(i).target_sz)) * DSKCFparameters.output_sigma_factor / DSKCFparameters.cell_size;
    scaleDSKCF_struct.yfs(i).yf = fft2(gaussian_shaped_labels( scaleDSKCF_struct.output_sigmas(i).output_sigma, floor( scaleDSKCF_struct.windows_sizes(i).window_sz / DSKCFparameters.cell_size)));
    
    %store pre-computed cosine window
    scaleDSKCF_struct.cos_windows(i).cos_window = hann(size(scaleDSKCF_struct.yfs(i).yf,1)) * hann(size(scaleDSKCF_struct.yfs(i).yf,2))';
    scaleDSKCF_struct.lens(i).len = scaleDSKCF_struct.target_sz(i).target_sz(1) * scaleDSKCF_struct.target_sz(i).target_sz(2);
    scaleDSKCF_struct.inds(i).ind = floor(linspace(1,scaleDSKCF_struct.lens(i).len, round(0.25 * scaleDSKCF_struct.lens(i).len)));
end

scaleDSKCF_struct.prevpos=pos;