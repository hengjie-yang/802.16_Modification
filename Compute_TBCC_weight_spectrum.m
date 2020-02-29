function weight_node = Compute_TBCC_weight_spectrum(constraint_length, code_generator, N)

% 
%   This function computes the weight spectrum of a given TBCC of length N
%
%   Inputs:
%       1) constraint_length: a scalar denoting the constraint length
%       2) code_generator: a matrix describing the CC generator
%       3) N: a scalar denoting the trellis depth (Note this is NOT the blocklength)
%
%   Outputs:
%       1)polys: a 2^v*1 symbolic vector, with the i-th entry denoting the
%       weight enumrating function of TBCC that starts at state i and end
%       at state i. i=1,2,...2^v, where v denotes the number of memory
%       elements.
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

% disp(T);

% Step 2: Compute the weight enumerating function for finite-length TBCC
disp('Step 2: compute the weight enumerating function for each starting state.');
B = eye(NumStates);
B = sym(B);

for iter = 1:N
    disp(['Current depths: ',num2str(iter)]);
    B = B*T;
    B = expand(B);
end
polys = diag(B);
% disp(polys);

% Step 3: Compute the overall weight enumerating function for finite-length
% TBCC
disp('Step 3: Compute the overall weight enumerating function.');

Poly = sum(polys);
disp(Poly);
weight_spectrum = coeffs(Poly,'All');
weight_spectrum = fliplr(weight_spectrum);
weight_spectrum = double(weight_spectrum);
weight_spectrum = weight_spectrum';

weight_node.weight_spectrum = weight_spectrum;
weight_node.overall_weight_function = Poly;
weight_node.weight_function_per_state = polys;


% Save results
code_string = '';
for iter = 1:size(code_generator,2)
    if iter < size(code_generator,2)
        code_string = [code_string, num2str(code_generator(iter)), '_'];
    else
        code_string = [code_string, num2str(code_generator(iter))];
    end
end

file_name = ['weight_node_TBCC_',code_string,'_N_',num2str(N),'.mat'];

save(file_name,'weight_node','-v7.3');




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
    

