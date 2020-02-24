% This script plots the RCU bound for k = 64, n = 2(k+32) = 192
%
%

set(0,'DefaultTextFontName','Times','DefaultTextFontSize',16,...
    'DefaultAxesFontName','Times','DefaultAxesFontSize',16,...
    'DefaultLineLineWidth',1,'DefaultLineMarkerSize',7.75);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');


SNR = -1:0.1:1.5;
k = 64;
m = 32;
tbcc_rate = 1/2;
N = (k+m)/tbcc_rate;
Rate = k/N;
Pe = zeros(1,size(SNR,2));

for iter = 1:size(SNR,2)
    snr = 10^(SNR(iter)/10);
    noise_var = 1/snr;
    [Pe(iter),~] = rcu(N, Rate, [0.5;0.5], [-1;1], noise_var);
end

semilogy(SNR,Pe,'-.','Color','black','LineWidth',2.0);
grid on
ylim([10^-6,1]);
xlabel('SNR (dB)','interpreter','latex');
ylabel('Frame Error Rate','interpreter','latex');
title('RCU bound for $k = 64$, $N = 192$, $R = 1/3$','interpreter','latex');

set(gcf,'Color','none');
export_fig RCU_bound_k_64_N_192.pdf -transparent

