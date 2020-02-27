function upper_bound = Search_upper_bound_d_CRC(code_generator, N, m)

% This function is to compute the upper bound of d_CRC; The upper bound is
% given by the distance at which the accumulative number of codewords first
% exceeds 2^(m-1).
%
%   Inputs:
%       1) code_generator: the matrix denoting the CC
%       2) N: the trellis depth
%       3) m: the CRC degree
%
%   Outputs:
%       1)upper bound: a scalar denoting the distance at which the
%       accumulative number of codewords first exceeds 2^(m-1).
%
%

%   Copyright 2020 Hengjie Yang

code_string = '';
for iter = 1:size(code_generator,2)
    if iter < size(code_generator,2)
        code_string = [code_string, num2str(code_generator(iter)), '_'];
    else
        code_string = [code_string, num2str(code_generator(iter))];
    end
end

file_name = ['weight_node_TBCC_',code_string,'_N_',num2str(N),'.mat'];
if ~exist(file_name, 'file')
    disp(['Error: the file ',file_name, ' does not exist!']);
    return
end

load(file_name, 'weight_node');
upper_bound = -1;
spectrum = weight_node.weight_spectrum;
maxDist = length(spectrum);
numCodes = 0;
threshold = 2^(m-1);
for dist = 1:maxDist
    if numCodes < threshold && numCodes + spectrum(dist) >= threshold
        upper_bound = 2*(dist - 1);
        break
    end
    numCodes = numCodes + spectrum(dist);
end


    
    



