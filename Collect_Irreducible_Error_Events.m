function IEE = Collect_Irreducible_Error_Events(constraint_length, code_generator, d_tilde, V, MaxLen)

% 
%   The collection algorithm collects all irreducible error events (IEEs)
%   of distance < d_tilde on the trellis of a rate-1/n tail-biting convolutional code (TBCC).
%   
%   Inputs:
%       1) constraint_length: a scalar indicating the constraint length of
%       TBCC
%       2) code_generator: a matrix specifying the generator of TBCC
%       3) d_tilde: the distance threshold
%       4) V: the predetermined ordering of states. The default ordering is
%       the natural ordering. The index starts from 1.
%
%   Outputs:
%       1) IEE is a structure with the following fields:
%           (1) state_ordering: a 1*2*v vector indicating the ordering of start
%           states.
%           (1) state_spectrum: a 2^v*1 vector that stores the aggregate of IEEs
%                   starting at each state specified by V.
%           (2) distance_spectrum: a d_tilde*1 vector that stores the
%                   aggregate of the IEEs with the given distance.
%           (2) list: a 2^v*1 cell, with the i-th cell representing the
%                   list of input sequence starting at state i.
%           (3) lengths: a 2^v*1 cell of vectors, with the i-th vector
%                   representing the list of lengths associated with the IEEs
%                   starting at state i.
%
%

%   Copyright 2020 Hengjie Yang



% Build the tail-biting trellis
trellis = poly2trellis(constraint_length, code_generator);
NumStates = trellis.numStates;
NumInputs = trellis.numInputSymbols; % for rate-1/n, NumInputs = 2

if nargin <= 3
    V = 1:NumStates; % state indices start from 1
end

if nargin <=4
    MaxLen = 200;
end




%Initialize some important parameters
IEE_list = cell(NumStates, 1);
IEE_lengths = cell(NumStates, 1);

inv_V = zeros(1, size(V,2));
for i = 1: NumStates
    inv_V(V(i)) = i;
end

for iter = 1:NumStates
        IEE_list{iter}=cell(d_tilde,1);
        IEE_lengths{iter}=cell(d_tilde, 1);
end


% Viterbi search at each starting state
for iter = 1:NumStates
    start_state = V(iter); % determine the start state
    disp(['Current progress: iter = ', num2str(iter), ' V(iter) = ',num2str(start_state)]);
    Column{1} = cell(NumStates, 1);
    Column{2} = cell(NumStates, 1);
    for depth = 0:MaxLen    % start a Viterbi search
        disp(['Current progress: iter = ', num2str(iter),...
            ' V(iter) = ',num2str(start_state),...
            ' Depth = ', num2str(depth)]);
        Column{mod(depth+1,2)+1}=cell(NumStates,1);
        if depth == 0
            for input = 1:NumInputs
                next_state = trellis.nextStates(start_state, input) + 1;
                if inv_V(next_state) >= iter % check if the state has been traversed in V
                    if isempty(Column{mod(depth+1, 2)+1}{next_state})
                        Column{mod(depth+1, 2)+1}{next_state} = cell(d_tilde, 1);%initialization
                    end
                    weight=sum(dec2bin(oct2dec(trellis.outputs(start_state, input)))-'0');
                    Column{mod(depth+1,2)+1}{next_state}{weight+1} = [Column{mod(depth+1,2)+1}{next_state}{weight+1};input-1]; 
                    % we use the assumption that NumInputs=2
                end
            end
        else
            for cur_state=1:NumStates
                if ~isempty(Column{mod(depth,2)+1}{cur_state})
                    if cur_state == start_state % meet the error event condition, stop extending
                        for dist = 1:d_tilde
                            if ~isempty(Column{mod(depth,2)+1}{cur_state}{dist})
                                [row_dim, col_dim] = size(IEE_list{start_state}{dist});
                                [row, col] = size(Column{mod(depth,2)+1}{cur_state}{dist});
                                if col_dim < col
                                    IEE_list{start_state}{dist} = [IEE_list{start_state}{dist}, Inf(row_dim, col-col_dim)];
                                else
                                    Column{mod(depth,2)+1}{cur_state}{dist} = [Column{mod(depth,2)+1}{cur_state}{dist}, Inf(row, col_dim-col)];
                                end
                                IEE_list{start_state}{dist} = [IEE_list{start_state}{dist}; Column{mod(depth,2)+1}{cur_state}{dist}];
                                IEE_lengths{start_state}{dist} = [IEE_lengths{start_state}{dist}; depth*ones(row, 1)];
                            end
                        end
                    else
                       for dist = 1:d_tilde
                           if ~isempty(Column{mod(depth,2)+1}{cur_state}{dist})
                               for input = 1:NumInputs
                                   next_state = trellis.nextStates(cur_state, input) + 1;
                                   if inv_V(next_state) >= iter
                                       if isempty(Column{mod(depth+1,2)+1}{next_state})
                                            Column{mod(depth+1,2)+1}{next_state} = cell(d_tilde, 1);
                                       end
                                       weight=sum(dec2bin(oct2dec(trellis.outputs(cur_state, input)))-'0');
                                       temp = Column{mod(depth,2)+1}{cur_state}{dist};
                                       temp = [temp, (input-1)*ones(size(temp,1), 1)];
                                       if dist + weight <= d_tilde
                                            Column{mod(depth+1,2)+1}{next_state}{dist+weight} = ...
                                                [Column{mod(depth+1,2)+1}{next_state}{dist+weight}; temp];
                                       end
                                   end
                               end
                           end
                       end
                    end
                end
            end
        end
    end
end

% Aggregate information
IEE.state_ordering = V;
IEE.state_spectrum = zeros(NumStates,1);
IEE.distance_spectrum = zeros(d_tilde,1);
for iter = 1:NumStates
    start_state = V(iter);
    for dist = 1:d_tilde
        IEE.state_spectrum(start_state) = IEE.state_spectrum(start_state) + length(IEE_lengths{start_state}{dist});
    end
end

for dist = 1:d_tilde
    for iter = 1:NumStates
        start_state = V(iter);
        IEE.distance_spectrum(dist) = IEE.distance_spectrum(dist) + length(IEE_lengths{start_state}{dist});
    end
end

IEE.list = IEE_list;
IEE.lengths = IEE_lengths;


code_string = '';
for iter = 1:size(code_generator,2)
    code_string = [code_string, num2str(code_generator(iter)), '_'];
end

save(['IEEs_TBCC_',code_string,'d_',num2str(d_tilde),'.mat'],...
    'IEE','-v7.3');


    
                                
                                
                                
                            
                    
            
            
        
        











