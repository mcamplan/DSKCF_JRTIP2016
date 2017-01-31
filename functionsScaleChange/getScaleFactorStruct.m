function scaleDSKCF_struct=getScaleFactorStruct(estimatedDepthMode,scaleDSKCF_struct) 
% GETSCALEFACTORSTRUCT.m select the target's scale
% 
%   GETSCALEFACTORSTRUCT this functions allows to select the new scales for
%   the DS-KCF tracker's model, according to the actual depth distribution
%   of the target. For more information about scale selection please see
%   [1]
%
%   INPUT:
%  -estimatedDepthMode depth of the target 
%  -scaleDSKCF_struct scale data structure (see INITDSKCFPARAM)
%  
%  OUTPUT
%  -scaleDSKCF_struct updated scale data structure
%
%
%  See also INITDSKCFPARAM
%
%  University of Bristol 
%  Massimo Camplani and Sion Hannuna
%  
%  massimo.camplani@bristol.ac.uk 
%  hannuna@compsci.bristol.ac.uk

scaleDSKCF_struct.updated = 0;

mode1 = estimatedDepthMode;
scaleDSKCF_struct.currDepth = mode1;

sf = scaleDSKCF_struct.InitialDepth / mode1;

% Check for significant scale difference to current scale
scaleOffset =  sf - scaleDSKCF_struct.scales(scaleDSKCF_struct.i);
if abs(scaleOffset) > scaleDSKCF_struct.minStep %% Need to change scale if possible
    if scaleOffset < 0 && scaleDSKCF_struct.i > 1% Getting smaller + check not smallest already
        diffs = scaleDSKCF_struct.scales(1:scaleDSKCF_struct.i) - sf;
        [a ind] = min(abs(diffs));
        if ind ~= scaleDSKCF_struct.i
            scaleDSKCF_struct.iPrev = scaleDSKCF_struct.i;
            scaleDSKCF_struct.i = ind;
            scaleDSKCF_struct.updated = 1;
        end
    elseif  scaleOffset > 0 && scaleDSKCF_struct.i < length(scaleDSKCF_struct.scales) % Getting bigger+ check not at biggest already
        diffs = scaleDSKCF_struct.scales(scaleDSKCF_struct.i:end) - sf;
        [a ind] = min(abs(diffs));
        ind = ind + scaleDSKCF_struct.i - 1;
        if ind ~= scaleDSKCF_struct.i
            scaleDSKCF_struct.iPrev = scaleDSKCF_struct.i;
            scaleDSKCF_struct.i = ind;
            scaleDSKCF_struct.updated = 1;
        end
    end  
end


