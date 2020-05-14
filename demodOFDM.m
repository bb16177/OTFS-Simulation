function [dataOut] = demodOFDM(dataIn,cpLen,ofdmSym)

%--------------------------------------------------------------------------
%
%              Demodulates OFDM symbols into serial data
%
%--------------------------------------------------------------------------
% Definition of input arguments 
%
% dataIn                         Input data vector
% cpLen                          Length of the cyclic prefix
% ofdmSym                        No. of ofdm symbols per subframe
%
%--------------------------------------------------------------------------
% Function returns: 
%
% dataOut                        Output time domain symbols
%
%--------------------------------------------------------------------------
%
% Author: Bradley Bates
% University of Bristol, UK
% email address: bb16177@bristol.ac.uk
% May 2020
%
% Code and algorithm originally from:
%
% Baher Mohammed (2020). OFDM signal generation, transmission and reception 
% (https://www.mathworks.com/matlabcentral/fileexchange/28368-ofdm-signal-generation-transmission-and-reception)
% MATLAB Central File Exchange. Retrieved May 10, 2020.
%
%--------------------------------------------------------------------------

% OFDM receiever reshapes serial data to parallel
parallelRx = reshape(dataIn, numel(dataIn)/ofdmSym, ofdmSym);
% Removing the cyclic Prefix
parallelRx(1:(cpLen), :) = [];

% Perform FFT
dataOut =  fft(parallelRx,[],2); 


end