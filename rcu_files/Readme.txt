The function signature is:

[Pe,rhop] = rcu(n, R, pX, X, sigma2)

n: Blocklength
R: spectral efficiency, i.e., log2(length(X))*R_c
pX: distribution on the constellation symbols
X: constellation symbols
sigma2: variance of AWGN noise

A typical call for your scenario to get the FER for rate 1/2 code with 
blocklength 128 bits on a biAWGN may look like:

Pe = rcu(128, 0.5, [0.5;0.5], [-1;1], 10^(-SNR/10));

If you need the other perspective, i.e., given a FER requirement (say 
1e-4), what is the necessary SNR, you may use this wrapper:

fzero(@(X) rcu(128, 0.5, [0.5;0.5], [-1;1], 10^(-X/10))-1e-4, [1 4]); % 
here, the interval [1,4] denotes a search range
