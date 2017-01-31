
 DS-KCF RGBD tracker code (version 1.0)

S. Hannuna, M. Camplani, J. Hall, M. Mirmehdi, D. Damen, T. Burghardt, A.
Paiement, L. Tao, DS-KCF: A ~real-time tracker for RGB-D data, Journal of
Real-Time Image Processing

doi:10.1007/s11554-016-0654-3

________________
This MATLAB code implements the DS-KCF tracker. 
This code is licensed under the BSD license

It is free for research use. If you find it useful, please acknowledge the paper
above with a reference. 


__________
Quickstart

1. Extract code somewhere. 

2. The tracker is ready to run on any of the 5 video sequences contained in
the data folder. 

DATA FOLDER CAN BE DOWNLOADED FROM
http://dx.doi.org/10.5523/bris.16vbnj3im1ygi1sh0yd0mt4lp0

or use the Rotational dataset in 
http://dx.doi.org/10.5523/bris.f8amc39n5f751ekgp2kbws1hb

3. Execute runDSKCF.m located in the testDS-KCFScripts folder 

4. Results will be generated in the DS-KCFresults folder


Documentation is located in the ds-kcfDoc folder. Other information and 
references to other toolboxes used can be found in Contents.m

Note: The tracker uses the 'fhog'/'gradientMex' functions from Piotr's Toolbox.
Some pre-compiled MEX files are provided for convenience. If they do not work for your
system, just get the toolbox from http://vision.ucsd.edu/~pdollar/toolbox/doc/index.html

Note: we tested this code in both windows and linux machine. Problems have 
been found with Matlab version 2015 while saving image results. Deactivate
image saving options if you are using that version (see runDSKCF.m) 
____________________________________

Copyright (c) 2016, Massimo Camplani, Sion Hannuna 
All rights reserved.
 
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
  1. Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
  2. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
 
  THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
  OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
  OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
  SUCH DAMAGE.
 
  