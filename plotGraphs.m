function [] = plotGraphs(berOFDM, berCOFDM, berOTFS, berCOTFS, M, numSC, EbNo)

%--------------------------------------------------------------------------
%
%           Plots and formats BER graphs of the input BER vectors
%
%--------------------------------------------------------------------------
% Input arguments: 
% 
% berOFDM                         OFDM ber vector
% berCOFDM                        coded-OFDM ber vector
% berOTFS                         OTFS ber vector
% berCOTFS                        coded-OFDM ber vector
% M                               Modulation order
% numSC                           no. subcarriers used in system
%
%--------------------------------------------------------------------------
% Function returns: 
% 
% void
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

% Set error rates of 0 to a very low value which can be plotted on graph
 berCOFDM(~berCOFDM)=1e-12;
 berCOTFS(~berCOTFS)=1e-12;

% Plotting 
figure

% Plot non-coded results
semilogy(EbNo,berOFDM(:,1),'-g');             %Plot simulated BER w/ OFDM
hold on;
semilogy(EbNo,berOTFS(:,1),'-r');             %Plot simulated BER w/ OTFS
% Plot coded results
semilogy(EbNo,berCOFDM(:,1),'--g');          %Plot simulated BER w/ C-OFDM
semilogy(EbNo,berCOTFS(:,1),'--r');             %Plot simulated BER w/ C-OTFS
    
% Formatting graph
axis([0 30 0.0005 0.5])
title(['BER for ', num2str(M), '-QAM in Wideband Rayleigh Fading with ', num2str(numSC),' Subcarriers']);
ylabel('Bit Error Rate');
xlabel('EbNo (dB)');
legend('OFDM', 'OTFS','C-OFDM','C-OTFS');
grid on;
hold off;


end