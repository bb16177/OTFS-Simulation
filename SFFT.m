function [outSig] = SFFT(inSig)

%--------------------------------------------------------------------------
%
%           Performs The Symplectic Fast Fourier Transform
%
%--------------------------------------------------------------------------
% Input arguments: 
% inSig                     Input N x M matrix to be transformed
%--------------------------------------------------------------------------
% Function returns: 
% outSig                    Output N x M matrix
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
outSig = sqrt(M/N) * ifft( fft(inSig, [], 1), [], 2);      % Apply transform

end