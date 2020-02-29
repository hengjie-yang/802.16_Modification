function storage = Estimate_memory_storage(constraint_length, code_generator, d_tilde, N)

% 
%   This function estimates the RAM storage required for "Reconstruct_TBP.m"
%
%   Inputs:
%       1) constraint_length: a scalar denoting the constraint length
%       2) code_generator: a matrix describing the CC generator
%       3) d_tilde: a scalar denoting the distance upper bound
%       3) N: a scalar denoting the trellis depth (Note this is NOT the blocklength)
%
%   Outputs:
%       1) storage: a scalar indicating the estimate of RAM storage in GB
%

%   Copyright 2020 Hengjie Yang


polys = [];
trellis = poly2trellis(constraint_length, code_generator);
NumStates = trellis.numStates;
T = zeros(NumStates, NumStates);
T = sym(T);
syms X;


% Step 1: compute the one-step transfer function
disp('Step 1: Compute the one-step transfer function.');
for cur_state = 1:NumStates
    for input = 1:2
        next_state = trellis.nextStates(cur_state,input) + 1;
        output_weight = sum(dec2bin(oct2dec(trellis.outputs(cur_state, input)))-'0');     
        output_symbol = convert_to_symbol(output_weight);
        T(cur_state, next_state) = output_symbol;
    end
end

% Step 2: Compute the weight enumerating function for all finite-length
% TBCCs
disp('Step 2: compute the weight enumerating function for each starting state.');
B = eye(NumStates);
B = sym(B);


total = 0;
for iter = 1:N
    disp(['Current depths: ',num2str(iter)]);
    B = B*T;
    B = expand(B);
    polys = diag(B);
    Poly = sum(polys); % weight enumrating function of length 'iter'
    weight_spectrum = coeffs(Poly,'All');
    weight_spectrum = fliplr(weight_spectrum);
    weight_spectrum = double(weight_spectrum);
    weight_spectrum = weight_spectrum';
    d_max = min(d_tilde, length(weight_spectrum));
    total = total + sum(weight_spectrum(1:d_max))*iter;
end

storage = total*(1e-9);
disp(['Estimate of memory: ', num2str(storage), ' GB']);




function chr = convert_to_symbol(weight)

chr = 0; %invalid string
syms X

if weight == 0
    chr = 1;
elseif weight == 1
    chr = X;
elseif weight == 2
    chr = X^2;
end





