function kf = linear_correlation(xf, yf)
%LINEAR_CORRELATION Linear Kernel at all shifts, i.e. correlation.
%   Computes the dot-product for all relative shifts between input images
%   X and Y, which must both be MxN. They must also be periodic (ie.,
%   pre-processed with a cosine window). The result is an MxN map of
%   responses.
%
%   Inputs and output are all in the Fourier domain.
%   This function has been inserted in the DS-KCF matlab library from the
%   KCF library released by
%
%   See also MAXRESPONSEDSKCF, MAXRESPONSEDEPTHWEIGHTDSKCF,
%   MODELUPDATEDSKCF, GAUSSIAN_CORRELATION, POLYNOMIAL_CORRELATION,  
%
%   Joao F. Henriques, 2014
%   http://www.isr.uc.pt/~henriques/
	
	%cross-correlation term in Fourier domain
	kf = sum(xf .* conj(yf), 3) / numel(xf);
    
end

