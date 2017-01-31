function [model_alphaf, model_alphaDf, model_xf, model_xDf]=modelUpdateDSKCF(...
    firstFrame,patch,patch_depth,features,cell_size,cos_window,kernel,yf,...
    lambda,model_alphaf, model_alphaDf, model_xf, model_xDf,scaleUpdate,interp_factor)
%MODELUPDATEDSKCF function for updating the DSKCF model %
%
%MODELUPDATEDSKCF.m is a function used for updating the DSKCF model according
%to the feature used. For more information about the model update see [1,2].
%Please note that  this function was partially built extending the KCF
%tracker code presented by Joao F. Henriques, in http://www.isr.uc.pt/~henriques/.
%
%
%  INPUT:
%  -patch patch of the color data
%  -patch_depth patch of the depth data
%  -features struct containing feature info.
%  -cell_size HOG parameter
%  -cos_window cosine window to smooth data in the Fourier domain
%  -kernel struct containing kernel information
%  -yf training labels
%  -lambda
%  - model_alphaf, model_alphaDf, model_xf, model_xDf are the models to be
%  updated (for depth and color if the two features are used indipendently)
%  -scaleUpdate flag set to one when the models templates need rescaling
%  -interp_factor interpolation factor for the model update
%
%  OUTPUT
%  - model_alphaf, model_alphaDf, model_xf, model_xDf are the updated
%  models (for depth and color if the two features are used indipendently)
%
% See also GAUSSIAN_CORRELATION, GET_FEATURES_DEPTH, UPDOWNSAMPLE_FOURIER
% POLYNOMIAL_CORRELATION, LINEAR_CORRELATION, SINGLEFRAMEDSKCF,
% SINGLEFRAMEDSKCF_OCCLUDER
%
%  [1] S. Hannuna, M. Camplani, J. Hall, M. Mirmehdi, D. Damen, T.
%  Burghardt, A. Paiement, L. Tao, DS-KCF: A real-time tracker for RGB-D
%  data, Journal of Real-Time Image Processing
%
%  [2] J. F. Henriques, R. Caseiro, P. Martins, and J. Batista. High-speed
%  tracking with kernelized correlation filters. Pattern Analysis and
%  Machine Intelligence, IEEE Transactions on, 2015.
%
%  University of Bristol
%  Massimo Camplani and Sion Hannuna
%
%  massimo.camplani@bristol.ac.uk
%  hannuna@compsci.bristol.ac.uk

%feature selection
if(features.hog_linear)
    
    %extract the features
    [xf xDf]=get_features_depth(patch,patch_depth, features, cell_size, cos_window);
    xf=fft2(xf);
    xDf=fft2(xDf);
    
    switch kernel.type
        case 'gaussian',
            kf = gaussian_correlation(xf, xf, kernel.sigma);
            kDf = gaussian_correlation(xDf, xDf, kernel.sigma);
        case 'polynomial',
            kf = polynomial_correlation(xf, xf, kernel.poly_a, kernel.poly_b);
            kDf = polynomial_correlation(xDf, xDf, kernel.poly_a, kernel.poly_b);
        case 'linear',
            kf = linear_correlation(xf, xf);
            kDf = linear_correlation(xDf, xDf);
    end
    
    alphaf = yf ./ (kf + lambda);
    alphaDf = yf ./ (kDf + lambda); %equation for fast training
    
    if (firstFrame),  %first frame, train with a single image
        model_alphaf = alphaf;
        model_alphaDf = alphaDf;
        model_xf = xf;
        model_xDf = xDf;
    else
        %subsequent frames, scale and interpolate model
        if(scaleUpdate)
            model_alphaf = updownsample_fourier( model_alphaf,size(alphaf,2),size(alphaf ,1));
            model_alphaDf = updownsample_fourier( model_alphaDf,size(alphaf,2),size(alphaf ,1));
            model_xf_ = xf;
            model_xDf_ = xDf;
            for i = 1:size(xf,3)
                model_xf_(:,:,i) = updownsample_fourier( model_xf(:,:,i),size(xf,2),size(xf ,1));
                model_xDf_(:,:,i) = updownsample_fourier( model_xDf(:,:,i),size(xDf,2),size(xDf ,1));
            end
            model_xf =  model_xf_;
            model_xDf =  model_xDf_;
            
            model_alphaf = (1 - interp_factor) * model_alphaf + interp_factor * alphaf;
            model_xf = (1 - interp_factor) * model_xf + interp_factor * xf;
            
            model_alphaDf = (1 - interp_factor) * model_alphaDf + interp_factor * alphaDf;
            model_xDf = (1 - interp_factor) * model_xDf + interp_factor * xDf;
        else
            %subsequent frames, interpolate model
            model_alphaf = (1 - interp_factor) * model_alphaf + interp_factor * alphaf;
            model_xf = (1 - interp_factor) * model_xf + interp_factor * xf;
            
            model_alphaDf = (1 - interp_factor) * model_alphaDf + interp_factor * alphaDf;
            model_xDf = (1 - interp_factor) * model_xDf + interp_factor * xDf;
        end
    end
    
else
    
    [xf ~]=get_features_depth(patch,patch_depth, features, cell_size, cos_window);
    xf=fft2(xf);
    %Kernel Ridge Regression, calculate alphas (in Fourier domain)
    switch kernel.type
        case 'gaussian',
            kf = gaussian_correlation(xf, xf, kernel.sigma);
        case 'polynomial',
            kf = polynomial_correlation(xf, xf, kernel.poly_a, kernel.poly_b);
        case 'linear',
            kf = linear_correlation(xf, xf);
    end
    alphaf = yf ./ (kf + lambda);   %equation for fast training
    
    if (firstFrame),  %first frame, train with a single image
        model_alphaf = alphaf;
        model_xf = xf;
    else
        
        %subsequent frames, interpolate model
        if(scaleUpdate)
            
            bigger_flag = 0;
            noChangeFlag= 0;
            if size(alphaf ,1) > size(model_alphaf, 1)
                %'bigger'
                bigger_flag = 0;
            elseif size(alphaf ,1) == size(model_alphaf, 1)
                %'should not be here'
                noChangeFlag=1;
            else
                %'smaller'
            end
            
            if bigger_flag%this is always zero in the final version....
                model_alphaf = updownsample_fourier_leave_gain( model_alphaf,size(alphaf,2),size(alphaf ,1), alphaf);
                model_xf_ = xf;
                for i = 1:size(xf,3)
                    model_xf_(:,:,i) = updownsample_fourier_leave_gain( model_xf(:,:,i),size(xf,2),size(xf ,1), xf(:,:,i));
                end
                model_xf =  model_xf_;
                
            elseif(noChangeFlag)
                aaaaa=0;
            else
                model_alphaf = updownsample_fourier( model_alphaf,size(alphaf,2),size(alphaf ,1));
                model_xf_ = xf;
                for i = 1:size(xf,3)
                    model_xf_(:,:,i) = updownsample_fourier( model_xf(:,:,i),size(xf,2),size(xf ,1));
                end
                model_xf =  model_xf_;
            end
            
            if bigger_flag
                % Don't do anything - already mixed
            else
                model_alphaf = (1 - interp_factor) * model_alphaf + interp_factor * alphaf;
                model_xf = (1 - interp_factor) * model_xf + interp_factor * xf;
            end
            
        else
            model_alphaf = (1 - interp_factor) * model_alphaf + interp_factor * alphaf;
            model_xf = (1 - interp_factor) * model_xf + interp_factor * xf;
        end
    end
    model_alphaDf=[];
    model_xDf=[];
end

end