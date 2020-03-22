function TBP_node = Reconstruct_TBPs(code_generator, d_tilde, N)

%
%   The function is a subfunction of "Search_DSO_CRC_by_Construction",
%   which only serves to the reconstruction of TBPs of length N. The
%   program is useful when designing a large degree CRC is required.
%
%   Inputs:
%       1) code_generator: a matrix specifying the generator of TBCC
%       2) d_tilde: the distance threshold
%       3) N: the trellis length
%
%   Outputs:
%       1) TBP_node: a struct including the following members:
%           (1) list: a d_tilde*1 cell indicating the list of length-N TBPs
%                   arranged in increasing distances. The true distance is
%                   the corresponding index minus 1.
%           (2) aggregate: a scalar indicating the number of length-N TBPs
%                   of distance less than 'd_tilde'.
%
%   Notes:
%       1) Must run "Collect_Irreducible_Error_Events.m" first if IEEs are
%       not generated before
%       2) Must run "Compute_TBCC_weight_spectrum.m" first if weight_node is
%       not generated before
%       3) The distance index equals true distance + 1, whereas the length
%       index equals true length.
% 

%   Copyright 2020 Hengjie Yang

tic

TBP_node  = {};
code_string = '';
for iter = 1:size(code_generator,2)
    code_string = [code_string, num2str(code_generator(iter)), '_'];
end

% write status messages in a .txt file
file_name = ['status_log_recon_TBPs_TBCC_',code_string,'d_',num2str(d_tilde),...
    '_N_',num2str(N),'.txt'];
StateFileID = fopen(file_name,'w');



file_name = ['IEEs_TBCC_',code_string,'d_',num2str(d_tilde),'.mat'];
if ~exist(file_name, 'file')
    msg = ['Error: the file ',file_name, ' does not exist!'];
    disp(msg);
    fprintf(StateFileID, '%s\n', msg);
    return
end
load(file_name, 'IEE');



file_name2 = ['weight_node_TBCC_',code_string,'N_',num2str(N),'.mat'];
if ~exist(file_name2, 'file')
    msg = ['Error: the file ',file_name2, ' does not exist!'];
    disp(msg);
    fprintf(StateFileID, '%s\n', msg);
    return
end
load(file_name2, 'weight_node');


full_spectrum = weight_node.weight_spectrum;
d_max = size(full_spectrum, 1);
if d_tilde > d_max
    msg = 'Error: (d_tilde-1) exceeds the maximum possible distance!';
    disp(msg);
    fprintf(StateFileID, '%s\n', msg);
    return
end

V = IEE.state_ordering;
NumStates = length(V);
% v = log2(NumStates);
Temp_TBPs = cell(d_tilde, 1); % used to find TBPs at each state



% State_spectrum = zeros(NumStates,1); % the state spectrum for all TBP(st_state)
Valid_TBPs = cell(d_tilde,1); % stores TBPs of length equal to N



% Step 1: Rebuild the tail-biting paths (TBPs) via IEEs
% Warning: the true distance = dist - 1 because we manually add 1
msg = 'Step 1: Reconstruct length-N TBPs using dynamic programming.';
disp(msg);
fprintf(StateFileID, '%s\n', msg);
for iter = 1:NumStates % find TBPs from every possible start state
    for dist = 1:d_tilde
        Temp_TBPs{dist} = cell(N+1, 1);
    end

    start_state = V(iter);
    List = IEE.list{start_state};
    
    %save memory
    for dist = 1:d_tilde
        List{dist} = int8(List{dist});
    end
    
    Lengths = IEE.lengths{start_state};
    msg = ['    Current start state: ', num2str(start_state),' number of IEEs: ',...
        num2str(IEE.state_spectrum(start_state))];
    disp(msg);
    fprintf(StateFileID, '%s\n', msg);
    for dist = 0:d_tilde-1 % true distance
        msg = ['        current distance: ',num2str(dist)];
        disp(msg);
        fprintf(StateFileID, '%s\n', msg);
        for len = 1:N %true length
            for weight = dist:-1:0 % true IEE weight
                for ii = 1:size(List{weight+1},1)
                    l = Lengths{weight+1}(ii);
                    if  weight == dist && l == len
                        Temp_TBPs{dist+1}{len+1} = int8(Temp_TBPs{dist+1}{len+1}); % save memory
                        Temp_TBPs{dist+1}{len+1} =[Temp_TBPs{dist+1}{len+1}; List{weight+1}(ii,1:l)];    
                    elseif l < len && ~isempty(Temp_TBPs{dist-weight+1}{len-l+1})
                        [row, ~] = size(Temp_TBPs{dist-weight+1}{len-l+1});
                        Added_bits = repmat(List{weight+1}(ii,1:l), row ,1);
                        New_TBPs = [Temp_TBPs{dist-weight+1}{len-l+1}, Added_bits];
                        Temp_TBPs{dist+1}{len+1} = int8(Temp_TBPs{dist+1}{len+1}); % save memory
                        Temp_TBPs{dist+1}{len+1} = [Temp_TBPs{dist+1}{len+1}; New_TBPs];
                    end
                end
            end
        end
    end

    % After building, we need to stack newly found TBPs to existing
    % TBPs


    for dist = 1:d_tilde
        if ~isempty(Temp_TBPs{dist}{N+1})
            Valid_TBPs{dist} = int8(Valid_TBPs{dist}); %save memory
            Valid_TBPs{dist} = [Valid_TBPs{dist}; Temp_TBPs{dist}{N+1}];
        end
    end
end



clearvars  Temp_TBPs







% Step 2: Build all valid TBPs through circular shift
msg = 'Step 2: Build remaining TBPs through cyclic shift.';
disp(msg);
fprintf(StateFileID, '%s\n', msg);
fclose(StateFileID);


parfor iter = 1:d_tilde
    file_name = ['status_log_recon_TBPs_TBCC_',code_string,'d_',num2str(d_tilde),...
    '_N_',num2str(N),'.txt'];
    StateFileID = fopen(file_name,'a');
    [row, ~] = size(Valid_TBPs{iter});
    
    % hash table was defined here.
    HashTable = containers.Map;
    for ii  = 1:row
        cur_seq = Valid_TBPs{iter}(ii,:);
        key_cur_seq = binary_to_hex(cur_seq);
        HashTable(key_cur_seq) = 1;
        
%         start_state = bin2dec(num2str(cur_seq(end-v+1:end)))+1;
%         State_spectrum(start_state) = State_spectrum(start_state)+1;
    end

    for ii = 1:row
        cur_seq = Valid_TBPs{iter}(ii,:);  
%         start_state = bin2dec(num2str(cur_seq(end-v+1:end)))+1;
        Extended_seq = [cur_seq, cur_seq]; 
        for shift = 1:N-1
            cyclic_seq = Extended_seq(1+shift:N+shift);
            key_cyclic_seq = binary_to_hex(cyclic_seq);
            if isequal(cyclic_seq, cur_seq) % termination condition for cyclic shift
                break
            end
            if ~isKey(HashTable, key_cyclic_seq)
                Valid_TBPs{iter} = [Valid_TBPs{iter};cyclic_seq]; % find a new TBP
                HashTable(key_cyclic_seq) = 1;
%                 State_spectrum(start_state) = State_spectrum(start_state)+1;
            end
        end
        msg = ['Step 2 Progress: ','dist = ',num2str(iter-1),...
            '   # of found TBPs: ',num2str(size(Valid_TBPs{iter},1)),...
            ' out of total: ',num2str(full_spectrum(iter)),...
            '    completed: ',num2str(size(Valid_TBPs{iter},1)/full_spectrum(iter)*100),'%'];
        disp(msg);
        fprintf(StateFileID, '%s\n', msg);
        
        if size(Valid_TBPs{iter},1) == full_spectrum(iter) % termination condition
            break
        end     
    end
end


TBP_node.list = Valid_TBPs;

% Miscellaneous: Compute total number of TBPs found
aggregate = 0;
for dist = 1:d_tilde
    aggregate = aggregate +size(Valid_TBPs{dist},1);
end

% TBP_node.state_spectrum = State_spectrum;

TBP_node.aggregate = aggregate;




% Save results
file_name = ['status_log_recon_TBPs_TBCC_',code_string,'d_',num2str(d_tilde),...
    '_N_',num2str(N),'.txt'];
StateFileID = fopen(file_name,'a');
msg = 'Completed and save results!';
disp(msg);
fprintf(StateFileID, '%s\n', msg);


file_name = ['TBP_node_TBCC_',code_string,'d_',num2str(d_tilde),'_N_',num2str(N),'.mat'];
save(file_name,'TBP_node','-v7.3');

timing = toc;
msg = ['Execution time: ',num2str(timing),'s'];
disp(msg);
fprintf(StateFileID, '%s\n', msg);
fclose(StateFileID);

end



function str = binary_to_hex(binary_vec)

str = '';
Len = length(binary_vec);
for i=Len:-4:1
    shift = min(3, i-1);
    four_tuple = binary_vec(i-shift:i);
    hex = dec2hex(bin2dec(num2str(four_tuple)));
    str = [str, hex];
end
str = fliplr(str);
end



