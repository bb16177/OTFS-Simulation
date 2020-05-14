function [outSig] = ISFFT(inSig)

%--------------------------------------------------------------------------
%
%           Performs Inverse Symplectic Fast Fourier Transform
%
%--------------------------------------------------------------------------
% Input arguments: 
% inSig                     Input N x M matrix to be transformed
%--------------------------------------------------------------------------
% Function returns: 
% outSig                    Output N x M matrix of doppler-Delay domain symbols
%--------------------------------------------------------------------------
%
% Author: Bradley Bates
% University of Bristol, UK
% email address: bb16177@bristol.ac.uk
% May 2020
%
% Copyright (c) 2020, Bradley Bates
%
%--------------------------------------------------------------------------

[N, M] = size(inSig);                                      % Calculate N and M
outSig = sqrt(N/M) * fft( ifft(inSig, [], 1), [], 2);      % Apply inverse transform

end