function [outSig] = multipathChannel(cpSize, delta_f, inSig, velocity)

%--------------------------------------------------------------------------
%
%               Generates and encodes random binary data
%
%--------------------------------------------------------------------------
% Input arguments: 
%
% cpSize                        Size of the cyclic prefix
% delta_f                       OFDM subcarrier spacing (Hz)
% inSig                         Input signal dimensions used to generate channel dimensions
% totalBits                     The approximate total bits to be simulated by the system
% velocity                      Channel mobility in km/hr
% 
%--------------------------------------------------------------------------
% Function returns: 
% 
% outSig                        Channel impulse response (serial vector)
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


% Create N x M channel matrix 
[N, M]= size(inSig');                                     % Size of inSig is used to create channel model
n = zeros(1,N);                                          % delay_Doppler rows (doppler)
m = zeros(1,M);                                          % delay_Doppler cols (delay)
H = transpose(n).*m;                                     % Create matrix

% Generate Channel Parameter
maxDelayspread = 0.5*((cpSize)/delta_f);  % Calculate max delay spread with 1/2 rule pf thumb
L = round(2*maxDelayspread * M*delta_f);  % Calculate necesary paths from bandwidth
step = maxDelayspread/L;                  % calculate difference between delays
pathDelays = (0:step:maxDelayspread);     % Discrete even delays of L-path channel
range shuffle;                            % Shuffle random no. generator
avgPathGains_dB = -(randi([3,7],[L,1]));  % Generate  random path gains in dB
avgPathGains = 10.^(0.1*avgPathGains_dB); % Convert to linear

% Calculate Max Doppler Shift
v = velocity*1e3/3600;                   % Mobile speed (m/s)
fc = 3.5e9;                              % Carrier frequency
fd = round(v*fc/physconst('lightspeed'));% Maximum Doppler shift to nearest Hz

% Generate doppler spreads w/ Jake's model
for l=0:L-1
    Vi(l+1) = fd* cos( (2*pi*l)/(L-1) );
end

% Initialize channel variables
T = 1/delta_f;                  % unextended OFDM symbol period
Ts = (1+cpSize)/delta_f;        % OFDM symbol period)
Ti = pathDelays;                % Path delays
hi = avgPathGains;              % Path gains

% Create matrix representation of channel
for m=1:M               % Loop along the rows
    for n=1:N           % Loop down the cols
        
        for x=1:L       % Loop to sum terms in the channel memory     
            % Define terms of model
            expTerm = (-2*1i*(pi)) * ((m+M/2).*delta_f.*Ti(x) - Vi(x).*(n).*Ts);
            hiPrime = hi(x)*(1 + 1i*(pi).*Vi(x).*T);
            % Produce channel impulse response model
            H(n, m) = H(n, m) + exp(expTerm) * hiPrime;
            
        end
    end

end


% % 3D plot of channnel model
% figure
% absH = 10*log10(abs(H));
% mesh(absH)
% surf(absH)
% colormap(jet)    % change color map
% shading interp    % interpolate colors across lines and faces
% hold on;
% title(['Time-Frequency Plot of The Fading Channel']);
% xlabel('Frequency (MHz)');
% xticks([8 142 274 408 540])
% xticklabels({'3.496','3.498','3.500','3.502','3.504'})
% ylabel('Time (ms)');
% yticks([0 14 28 42 56 70 84 98 112 126 140])
% yticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
% zlabel('Power (dB');
% hold off;


% Convert to serial
outSig = reshape(H',[n*m,1]);
end
