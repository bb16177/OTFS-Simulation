function [dataOut] = modOFDM(dataIn,numSC,cpLen,ofdmSym)

%--------------------------------------------------------------------------
%
%               Modulates random input data into OFDM symbols
%
%--------------------------------------------------------------------------
% Input Arguments: 
% dataIn                         Input data vector
% numSC                          The number of subcarriers used
% cpLen                          Length of the cyclic prefix
% ofdmSym                        No. of ofdm symbols per subframe
%--------------------------------------------------------------------------
% Function returns: 
% dataOut                        Output OFDM symbols
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

% Calculate variables
cyclicPrefix_start  = numSC - cpLen;

% Perform IFFT
ifftSubcarrier = ifft(dataIn,[],2); 

%Finding cyclic prefix for each subcarrier
for i=1:cpLen
    for j=1:ofdmSym                   
        cyclicPrefix_data(i,j) = ifftSubcarrier(i+cyclicPrefix_start,j);
    end
end

% Add cyclic prefix to the data
appendedCP = vertcat(cyclicPrefix_data, ifftSubcarrier);

% Convert to serial
dataOut = reshape(appendedCP,[numel(appendedCP),1]);



end
