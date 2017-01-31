function noiseStd=calculateNoiseVar(depthValue,  noiseModelP0, noiseModelP1,noiseModelP2)
% CALCULATENOISEVAR.m Calculates the variance of the Kinect Noise 
% 
%   CALCULATENOISEVAR is a function for calculating the depth noise for a
%   given distance according to the quadratic noise model presented in [1]
%
%   INPUT: 
%  -depthValue depth measurement
%  -noiseModelP0, noiseModelP1,noiseModelP2 noise model parameters
%
%   OUTPUT
%  -noiseSTD estimated noise standard deviation 
%  
%  See also FASTDEPTHSEGMENTATIONDSKCF_NOISEMODEL,
%  CHECKOCCLUSIONSDSKCF_NOISEMODEL, CHECKOCCLUSIONSDSKCF_SECONDPLANE
%
%  [1] M. Camplani, T. Mantecon, and L. Salgado. Depth-color fusion
%  strategy for 3-D scene modeling with Kinect. Cybernetics, IEEE
%  Transactions on, 43(6):1560–1571, 2013
%
%
%  University of Bristol 
%  Massimo Camplani and Sion Hannuna
%  
%  massimo.camplani@bristol.ac.uk 
%  hannuna@compsci.bristol.ac.uk

noiseStd=noiseModelP0+ noiseModelP1*depthValue+ noiseModelP2*depthValue*depthValue;
