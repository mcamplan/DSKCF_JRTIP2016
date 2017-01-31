function [L,Cnew,LUT,H,I,LUTCC]=fastDepthSegmentationDSKCF_noiseModel(im,c,nanMatrix,minimumError,Cinit, findPeak,targetDepth,targetSTD,noiseModelVector)
% FASTDEPTHSEGMENTATIONDSKCF_NOISEMODEL.m segments depth data
% 
%   FASTDEPTHSEGMENTATIONDSKCF_NOISEMODEL function applies the fast depth
%   segmentation algorithm described in [1]. The segmentation is composed
%   by two different stage a fast version of the Kmeans applied to the
%   depth data, plus a connected component analysis to refine clusters in
%   the image plane. This function was implemented by starting from the
%   "Fast segmentation of N-dimensional grayscale images" presented by
%   Anton Semechko and shared in the Matlab Central at this link under BSD
%   licence
%   http://www.mathworks.com/matlabcentral/fileexchange/41967-fast-segmentation-of-n-dimensional-grayscale-images
%
%   INPUT: 
%   - im   depth image coded in 16bits, each pixel contains mm data. 
%   - c    positive interger greater than 1 specifying the number of
%           clusters. c=2 is the default setting. Alternatively, c
%           initialized by considering peaks in the depth distribution.
%   -nanMatrix  binary mask containing flags for missing depth pixels 
%   -minimumError convergence criteria for the Kmeans algorithm
%   -histStep histogram bin used to compose depth histogram
%   -Cinit initial Kmeans seeds. Set this values to -1 to not initialize
%   externally the starting seeds
%   -findPeak boolean flag to initialize the Kmeans seed with the peaks in
%   the depth distribution 
%   -targetDepth,targetSTD  depth mean value and standard deviation of the
%   object in the previous frame. THey are used to calculate the noise
%   component
%
%   OUTPUT
%   - L    label image of the same size as the input image. For example,
%           L==i represents the region associated with prototype C(i),
%           where i=[1,k] (k = number of clusters).
%   - Cnew    1-by-k array of cluster centroids.
%   - LUT  L-by-1 array that specifies the intensity-class relations,
%           where L is the dynamic intensity range of the input image. 
%           Specifically, LUT(1) corresponds to class assigned to 
%           min(im(:)) and LUT(L) corresponds to the class assigned to
%           max(im(:)). 
%   -H histogram's bin height
%   -I histogram's bin centers
%   -LUTCC the same as LUT but considering also the connected component
%   anaylsis
%
%  See also LUT2LABELNANSUPPORT, LUT2LABELNANSUPPORTCC, CALCULATENOISEVAR,
%  CHECKOCCLUSIONSDSKCF_NOISEMODEL, SINGLEFRAMEDSKCF,
%  CHECKOCCLUSIONSDSKCF_SECONDPLANE
%
% [1] S. Hannuna, M. Camplani, J. Hall, M. Mirmehdi, D. Damen, T.
% Burghardt, A.Paiement, L. Tao, DS-KCF: A ~real-time tracker for RGB-D
% data, Journal of Real-Time Image Processing
%
%
%  University of Bristol 
%  Massimo Camplani and Sion Hannuna
%  
%  massimo.camplani@bristol.ac.uk 
%  hannuna@compsci.bristol.ac.uk

% Default input arguments
if nargin<2 || isempty(c), c=2; end

% Basic error checking
if nargin<1 || isempty(im)
    error('Insufficient number of input arguments')
end
msg='Revise variable used to specify class centroids. See function documentaion for more info.';
if ~isnumeric(c) || ~isvector(c)
    error(msg)
end
if numel(c)==1 && (~isnumeric(c) || round(c)~=c || c<2)
    error(msg)
end

% Check image format
if isempty(strfind(class(im),'int'))
    error('Input image must be specified in integer format (e.g. uint8, int16)')
end
if sum(isnan(im(:)))~=0 || sum(isinf(im(:)))~=0
    error('Input image contains NaNs or Inf values. Remove them and try again.')
end


histStep=max(2.5*calculateNoiseVar(targetDepth,noiseModelVector(1),noiseModelVector(2),noiseModelVector(3)),targetSTD);

newPointSet=im(~nanMatrix);
Imin=double(min(newPointSet));
Imax=double(max(newPointSet));
I=(Imin:histStep:Imax)';
if(isempty(I))
    L=[];
    Cnew=[];
    LUT=[];
    H=[];
    I=[];
    LUTCC=[];
    return
end

if(I(end)~=Imax || length(I)==1)
    I(end+1)=Imax+histStep;
end
I=I(:);
% Compute intensity histogram
H=hist(double(newPointSet),I);
H=H(:);
maxValue=max(H);

minPeakDistParam=3;
if(length(I)<50)
   minPeakDistParam=2;
end

[peakDepth,posPeak]=findpeaks([0; H ;0],'MINPEAKDISTANCE',minPeakDistParam,'MINPEAKHEIGHT',0.005*maxValue);

% Initialize cluster centroids
if numel(c)>1
    C=c;
    c=numel(c);
else
    dI=(Imax-Imin)/c;
    if(isempty(Cinit))
        C=Imin+dI/2:dI:Imax;
    else
        C=Cinit;
    end
end

if(findPeak)
    if(length(C)==length(posPeak) && C(1)==-1);
        C=I(posPeak-1);
    elseif (length(C)~=length(posPeak))
        c=length(posPeak);
        C=I(posPeak-1);
    end
    C=C';
end



% Update cluster centroids
IH=I.*H; dC=Inf;

C0=C;
Citer=C;
while dC>minimumError
    
    Citer=C;
    
    % Distance to the centroids
    D=abs(bsxfun(@minus,I,C));
    
    % Classify by proximity
    [Dmin,LUT]=min(D,[],2); %#ok<*ASGLU>
    for j=1:c
        C(j)=sum(IH(LUT==j))/sum(H(LUT==j));
        if(isnan(C(j)))
            C(j)=Citer(j);
        end
    end
      
    % Change in centroids 
    dC=max(abs(C-Citer));
    
end

%%correct singleton peaks
if(length(C)==1)
    C=I(posPeak-1);
end

[L,Cnew,LUTCC]=LUT2labelNanSupportCC(im,LUT,nanMatrix,histStep,C);

