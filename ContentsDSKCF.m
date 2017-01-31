% DS-KCF RGBD tracker code Info and Contents
%
% 
% S. Hannuna, M. Camplani, J. Hall, M. Mirmehdi, D. Damen, T. Burghardt, A.
% Paiement, L. Tao, DS-KCF: A ~real-time tracker for RGB-D data, Journal of
% Real-Time Image Processing
%
% DS-KCF main functions.
%   wrapperDSKCF   - DS-KCF tracker wrapper function.
%   singleFrameDSKCF - DS-KCF tracker core function 
%   modelUpdateDSKCF - DS-KCF tracker model update
%   occludingObjectSegDSKCF - occluding object segmentation
%
% DS-KCF test script
%   testDS-KCFScripts\runDSKCF.m    
%
% External toolboxes used.
%   KCF [1] matlab library presented by Joao F.Henriques, in
%   http://www.isr.uc.pt/~henriques/
%
%   Piotr's Toolbox
%   http://vision.ucsd.edu/~pdollar/toolbox/doc/index.html
%
%   RGBD tracker [2] and Benchmark
%   http://tracking.cs.princeton.edu/code.html
%
%  [1] J. F. Henriques, R. Caseiro, P. Martins, and J. Batista. High-speed
%  tracking with kernelized correlation filters. Pattern Analysis and
%  Machine Intelligence, IEEE Transactions on, 2015.
%
%  [2] Shuran Song and Jianxiong Xiao. Tracking Revisited using RGBD
%  Camera: Baseline and Benchmark. 2013.
%
% License: This MATLAB code implements the DS-KCF tracker. 
% This code is licensed under the BSD license
%
%  University of Bristol
%  Massimo Camplani and Sion Hannuna
%  
%  massimo.camplani@bristol.ac.uk
%  hannuna@compsci.bristol.ac.uk
%  Copyright (c) 2016, Massimo Camplani, Sion Hannuna
%  All rights reserved.
%  
%   Redistribution and use in source and binary forms, with or without
%   modification, are permitted provided that the following conditions are
%   met: 
%   1. Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%   2. Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in the
%      documentation and/or other materials provided with the distribution.
%  
%   THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
%   ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
%   PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS
%   BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
%   BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
%   WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
%   OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
%   ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.