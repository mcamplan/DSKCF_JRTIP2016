function [L,Cnew,LUTCC]=LUT2labelNanSupportCC(im,LUT,nanMatrix,histStep,C)
% LUT2LABELNANSUPPORTCC.m is a function for assign clustering label to the segmented depth image 
% 
%   LUT2LABELNANSUPPORTCC function creates the clusters'label for all the
%   pixels of the segmented depth image. A cluster refinement is performed
%   by performing a connected component anaylisis on the image plane. This
%   function extends the "Fast segmentation of N-dimensional grayscale
%   images" presented by Anton Semechko and shared in the Matlab Central at
%   this link under BSD licence
%   http://www.mathworks.com/matlabcentral/fileexchange/41967-fast-segmentation-of-n-dimensional-grayscale-images
%
%   INPUT: 
%   - im   depth image coded in 16bits, each pixel contains mm data. 
%   -LUT is the lookuptable containing cluster label and corresponding
%   depth value of the histogram
%   -nanMatrix  binary mask containing that marks missing depth pixels 
%   -histStep histogram bin used to compose depth histogram
%   -C centroids of the depth clusters
%
%   OUTPUT
%   - L    label image of the same size as the input image. For example,
%           L==i represents the region associated with prototype C(i),
%           where i=[1,k] (k = number of clusters).
%   -Cnew new centroids of the depth clusters obtained after the connected
%   component analysis
%
%   -LUTCC new LUT obtained after the connected component analysis
%
% See also FASTDEPTHSEGMENTATIONDSKCF_NOISEMODEL
%
%  University of Bristol 
%  Massimo Camplani and Sion Hannuna
%  
%  massimo.camplani@bristol.ac.uk 
%  hannuna@compsci.bristol.ac.uk

Cnew=[];
% Intensity range
newPointSet=im(~nanMatrix);
Imin=double(min(newPointSet));
Imax=double(max(newPointSet));
I=(Imin:histStep:Imax)';

offsetHist=histStep/2;

LUTCC=[];
%I(end)=Imax;
if(I(end)~=Imax|| length(I)==1)
    I(end+1)=Imax+histStep;
end


% Create label image
L=zeros(size(im),'uint8');
detectedRegion=0;
startingIndex=1;
for k=1:max(LUT)
   
    % Intensity range for k-th class
    i=find(LUT==k);
    if(isempty(i)==false)
        i1=i(1);
        if numel(i)>1
            i2=i(end);
        else
            i2=i1;
        end
        
        % Map the intensities in the range [I(i1),I(i2)] to class k
        bw=im>=I(i1)-offsetHist & im<I(i2)+offsetHist;
        [Ltemp num]=bwlabel(bw);
        Ltemp=Ltemp+detectedRegion;
        detectedRegion=detectedRegion+num;
        L(bw)=Ltemp(bw);
        Cnew(startingIndex:(startingIndex+num-1))=C(k);
        LUTCC(startingIndex:(startingIndex+num-1))=k;
        startingIndex=startingIndex+num;
    end
end

