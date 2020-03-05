x=[1024, 4095, 16341, 63179, 221547, 655426]; % number of TBPs
y=[9, 11.15, 13, 17.7, 42.9, 386.15]; % execution time
log_y = log10(y);
% estimated_y = 2.2e-6*x+1.15;

estimated_y = 10^1.15*10.^(2.2e-6*x);


plot(x,y,'-*');hold on
plot(x,estimated_y,'-o');
grid on