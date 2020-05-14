function [dataIn, dataBits_in, codedData_in, packetSize, numPackets, numCB] = dataGen(k,numDC,ofdmSym,totalBits,codeRate, ldpcEncoder)

%--------------------------------------------------------------------------
%
%               Generates and encodes random binary data
%
%--------------------------------------------------------------------------
% Input arguments: 
%
% k                             Bits/Symbol
% numDC                         Number of data subcarriers
% ofdmSym                       Totol ofdm symbols per subframe
% totalBits                     The approximate total bits to be simulated by the system
% codeRate                      LDPC code rate
% ldpcEncoder                   LDPC encode system object
% 
%--------------------------------------------------------------------------
% Function returns: 
% 
% dataIn                        The input encoded data to the modulator
% dataBits_in                   The binary generated initial data before coding
% codedData_in                  The binary codewords
% packetSize                    No. of subframes / "packet"
% numPackets                    Number of packets required to satisfy totalBits
% numCB                         No. of code blocks / "subframe"
%
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

% Calculate subframe size
% Initialise information about frames to be transmitted
packetSize = 1;              % No. of subframes / "packet"
numCB = 1;                   % No. of code blocks / "subframe"
noCodedbits = 64800;         % Codeword length

% Calculate exact no. of frames and bits required for simulation
subframeSize = [k*numDC*ofdmSym 1];               % Calculate size of a subframe
maxSubframes = ceil(totalBits./subframeSize(1));  % Total no. of subframes to be transmitted

% Determine number of code block and number of subframes/packet
if subframeSize(1) == noCodedbits               % If same size do nothing
    numCB = 1;
    packetSize = 1;
    
elseif subframeSize(1) > noCodedbits            % Match for when subframe > codeword
    [numCB, packetSize] = rat(subframeSize(1)./ noCodedbits,1e-1);
    
elseif subframeSize(1) < noCodedbits            % Match for when subframe < codeword
    [packetSize, numCB] = rat(noCodedbits./ subframeSize(1),1e-1);
end   

% Ensure theres always enough data bit capacity
while numCB*noCodedbits >= subframeSize(1)*packetSize
    packetSize = packetSize + 1;
    % Divide both by 2 if possible
    if (rem(numCB,2) == 0) && (rem(packetSize,2) == 0) && numCB*noCodedbits <= subframeSize(1)*packetSize
        packetSize = packetSize./2;
        numCB = numCB./2;
    end
end 

% Calculate the pad bits required to round up to the nearest whole subframe
padBits = zeros((subframeSize(1)*packetSize - numCB*noCodedbits),1);% Calculate the number of pad bits
numPackets = round(maxSubframes./packetSize);                       % Total no. of "packets" to be transmitted
numPackets(~numPackets)=1;                                          % If 0 packets set to 1

% Generate Random Input Data
range shuffle;                               % Shuffle random no. generator
codedData_in = [];                           % Initialse arrays
dataBits_in = [];
for q = 1:numCB
    dataBits = randi([0,1],[noCodedbits*codeRate,1]);      % Generate binary datafor each frame of the burst
    dataBits_in = [dataBits_in;dataBits];                  % Concatenate binary data  
    codedData_in = [codedData_in; ldpcEncoder(dataBits)];  % Generate LDPC codeblocks
end
paddedData_in = [codedData_in; padBits];             % Append pad bits to coded bits
dataIn = randintrlv(paddedData_in,4831);             % Interleave paddedData randdeintrlv


% Print info about data output
% fprintf(['\n',num2str(numCB),' code block(s)\n', num2str(packetSize), ' subframes/packet \n'])
% fprintf([ num2str(numel(padBits)),' pad bits out of ', num2str((subframeSize(1)*packetSize)),' data bits = ', num2str(100*(numel(padBits)/(subframeSize(1)*packetSize))),' percent pad bits!\n'])
% fprintf(['Total packets: ', num2str(numPackets),'\n'])
% fprintf(['Total bits:    ', num2str(numPackets*packetSize*(subframeSize(1))),'\n'])

end