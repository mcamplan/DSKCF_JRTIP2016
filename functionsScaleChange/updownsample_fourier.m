function out_m = updownsample_fourier( in_m,out_x_sz,out_y_sz)
% UPDOWNSAMPLE_FOURIER.m upsamples DS-KCF models in the Fourier Domain
%
%   UPDOWNSAMPLE_FOURIER is a function for upsampling DS-KCF model in the
%   Fourier Domain. This operation is fundamental for the DS-KCF tracker as
%   described in [1]. This function was implemented by using the function
%   "updownsample" presented by Ohad Gal and shared in the Matlab Central
%   at this link under BSD licence
%   http://uk.mathworks.com/matlabcentral/fileexchange/4658-updownsample
%
% input:    in_m                - input matrix for up/down sampling. can be in
%                                 space domain OR in fourier domain, in such
%                                 case, needs to be in matlab format !!!
%                                 (matlab-format = the save as given from fft/fft2)
%           out_x_sz,out_y_sz   - desired number of pixels in the output image
%           is_fourier_flag     - 1: the input is given in the fourier domain
%                                 0: the input is given in the space domain
%                                    (we need to use fft2 to convert to fourier domain)
%           is_real_flag        - 0: the input is a complex matrix -> don't use
%                                    abs() at the output, perform complex
%                                    up/down sampling
%                                 1: the input is real BUT has negative values ->
%                                    use real() at the output
%                                 2: the input is real and positive -> using
%                                    abs() at the output
%
% output:   out_m               - up/down sampled image
%
% See also MODELUPDATEDSKCF
%
%  [1] S. Hannuna, M. Camplani, J. Hall, M. Mirmehdi, D. Damen, T.
%  Burghardt, A. Paiement, L. Tao, DS-KCF: A real-time tracker for RGB-D
%  data, Journal of Real-Time Image Processing
%
%  University of Bristol
%  Massimo Camplani and Sion Hannuna
%
%  massimo.camplani@bristol.ac.uk
%  hannuna@compsci.bristol.ac.uk

% ==============================================
% get input image size, and calculate the gain
% ==============================================
[in_y_sz,in_x_sz] = size( in_m );
gain_x = out_x_sz/in_x_sz;
gain_y = out_y_sz/in_y_sz;


% build grid vectors for the up/down sampling
% ============================================
% if the input is even & output is odd-> use floor for all
% if the output is even & input is odd -> use ceil for all
% other cases - don't care
% for downsampling -> the opposite
if (~mod( in_x_sz,2 ) & (out_x_sz>in_x_sz)) | (mod( in_x_sz,2 ) & (out_x_sz<in_x_sz))
    x_output_space  = max(floor((out_x_sz-in_x_sz)/2),0) + [1:min(in_x_sz,out_x_sz)];
    x_input_space   = max(floor((in_x_sz-out_x_sz)/2),0) + [1:min(in_x_sz,out_x_sz)];
else
    x_output_space  = max(ceil((out_x_sz-in_x_sz)/2),0) + [1:min(in_x_sz,out_x_sz)];
    x_input_space   = max(ceil((in_x_sz-out_x_sz)/2),0) + [1:min(in_x_sz,out_x_sz)];
end
if (~mod( in_y_sz,2 ) & (out_y_sz>in_y_sz)) | (mod( in_y_sz,2 ) & (out_y_sz<in_y_sz))
    y_output_space  = max(floor((out_y_sz-in_y_sz)/2),0) + [1:min(in_y_sz,out_y_sz)];
    y_input_space   = max(floor((in_y_sz-out_y_sz)/2),0) + [1:min(in_y_sz,out_y_sz)];
else
    y_output_space  = max(ceil((out_y_sz-in_y_sz)/2),0) + [1:min(in_y_sz,out_y_sz)];
    y_input_space   = max(ceil((in_y_sz-out_y_sz)/2),0) + [1:min(in_y_sz,out_y_sz)];
end

% perform the up/down sampling
padded_out_m    = zeros( out_y_sz,out_x_sz );
in_m            = fftshift(in_m);
padded_out_m( y_output_space,x_output_space ) = in_m(y_input_space,x_input_space);
%out_m           = (gain_x*gain_y)*ifft2(ifftshift(padded_out_m));
out_m           = gain_x*gain_y*ifftshift(padded_out_m);



