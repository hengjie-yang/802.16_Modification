function Extract_network_from_trellis(constraint_length, code_generator)

%
%   This function is to automatically generate the network file from the
%   trellis.
%
%   Inputs:
%       1) constraint_length: a scalar denoting the constraint length
%       2) code_generator: a matrix describing the convolutional encoder
%
%   Outputs:
%       1) Network: a matrix specified by the mason file
%

%   Copyright 2020 Hengjie Yang

trellis = poly2trellis(constraint_length, code_generator);
NumStates = trellis.numStates;
Network = zeros(NumStates,2);
output_symbols = [];

id = 1;
for cur_state = 1:NumStates
    for input = 1:2
        next_state = trellis.nextStates(cur_state,input);
        output_weight = sum(dec2bin(oct2dec(trellis.outputs(cur_state, input)))-'0');     
        output_symbol = convert_to_symbol(output_weight);
        output_symbols = [output_symbols; output_symbol];
        Network(id,:) = [cur_state-1, next_state];
        id = id + 1;
    end
end

Indices = 1:2*NumStates;
Indices = Indices';


code_string = '';
for iter = 1:size(code_generator,2)
    if iter < size(code_generator,2)
        code_string = [code_string, num2str(code_generator(iter)), '_'];
    else
        code_string = [code_string, num2str(code_generator(iter))];
    end
end

filename = ['CC_',code_string,'.txt'];
fileID = fopen(filename, 'w');
data = [Indices, Network, output_symbols];
data = data';
fprintf(fileID, '%s\t %s\t %s\t %s\n', data);
fclose(fileID);


function chr = convert_to_symbol(weight)

chr = "-1"; %invalid string

if weight == 0
    chr = "(1)";
elseif weight == 1
    chr = 'X';
elseif weight == 2
    chr = 'X^2';
end
    


