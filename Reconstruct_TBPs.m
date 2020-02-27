function TBP_node = Reconstruct_TBPs(code_generator, d_tilde, N)

%
%   The function is a subfunction of "Search_DSO_CRC_by_Construction",
%   which only serves to the reconstruction of TBPs of length N. The
%   program is useful when designing a large degree CRC is required.
%
%   Inputs:
%       1) code_generator: a matrix specifying the generator of TBCC
%       3) d_tilde: the distance threshold
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
%       1) Must run "Collect_Irreducible_Error_Events" first if IEEs are
%       not generated before
%       2) The distance index equals true distance + 1, whereas the length
%       index equals true length.
% 



code_string = '';
for iter = 1:size(code_generator,2)
    code_string = [code_string, num2str(code_generator(iter)), '_'];
end

file_name = ['IEEs_CC_',code_string,'dtilde_',num2str(d_tilde),'.mat'];
if ~exist(file_name, 'file')
    disp(['Error: the file ',file_name, ' does not exist!']);
    return
end

load(file_name, 'IEE');


V = IEE.state_ordering;
NumStates = length(V);
Temp_TBPs = cell(d_tilde, 1); % used to find TBPs at each state


TBPs = cell(d_tilde, 1); % official list
for dist = 1:d_tilde
    TBPs{dist} = cell(N+1,1);
end

State_spectrum = zeros(NumStates,1);
Valid_TBPs = cell(d_tilde,1); % stores TBPs of length equal to N



% Step 1: Rebuild the tail-biting paths (TBPs) via IEEs
% Warning: the true distance = dist - 1 because we manually add 1
disp('Step 1: Reconstruct length-N TBPs using dynamic programming.');
for iter = 1:NumStates % find TBPs from every possible start state
    for dist = 1:d_tilde
        Temp_TBPs{dist} = cell(N+1, 1);
    end

    start_state = V(iter);
    List = IEE.list{start_state};
    Lengths = IEE.lengths{start_state};
    disp(['    Current start state: ', num2str(start_state),' number of IEEs: ',...
        num2str(IEE.state_spectrum(start_state))]);

    for dist = 0:d_tilde-1 % true distance
        for len = 1:N %true length
            for weight = dist:-1:0 % true IEE weight
                for ii = 1:size(List{weight+1},1)
                    l = Lengths{weight+1}(ii);
                    if  weight == dist && l == len
                        Temp_TBPs{dist+1}{len+1} =[Temp_TBPs{dist+1}{len+1}; List{weight+1}(ii,1:l)];    
                    elseif l < len && ~isempty(Temp_TBPs{dist-weight+1}{len-l+1})
                        [row, ~] = size(Temp_TBPs{dist-weight+1}{len-l+1});
                        Added_bits = repmat(List{weight+1}(ii,1:l), row ,1);
                        New_TBPs = [Temp_TBPs{dist-weight+1}{len-l+1}, Added_bits];
                        Temp_TBPs{dist+1}{len+1} = [Temp_TBPs{dist+1}{len+1}; New_TBPs];
                    end
                end
            end
        end
    end

    % After building, we need to stack newly found TBPs to existing
    % TBPs
    for dist = 1:d_tilde
        for len = 1:N+1
            if ~isempty(Temp_TBPs{dist}{len})
                TBPs{dist}{len} = [TBPs{dist}{len}; Temp_TBPs{dist}{len}];
            end
        end
    end
end



for dist=1:d_tilde
    if ~isempty(TBPs{dist}{N+1})
        Valid_TBPs{dist} = TBPs{dist}{N+1};
    end
end


clearvars TBPs Temp_TBPs


% % Verify if there are repetitive TBPs after building
% HashNumber = 2^N+1; % This is the maximum number of cyclic shift
% HashTable = zeros(HashNumber, 1);
% for iter = 1:d_tilde
%     for ii = 1:size(Valid_TBPs{iter},1)
%         cur_seq = Valid_TBPs{iter}(ii,:);
%         h = ComputeHash(cur_seq, HashNumber);
%         HashTable(h) = HashTable(h) + 1;
%     end
% end


% Step 2: Build all valid TBPs through circular shift
disp('Step 2: Build remaining TBPs through cyclic shift.');
parfor iter = 1:d_tilde
    disp(['    Current distance: ',num2str(iter-1)]);
    [row, ~] = size(Valid_TBPs{iter});
    % hash table was defined here.

    HashTable = containers.Map;
    for ii  = 1:size(Valid_TBPs{iter},1)
        cur_seq = Valid_TBPs{iter}(ii,:);
        key_cur_seq = dec2bin(bi2de(cur_seq),N);
        HashTable(key_cur_seq) = 1;
    end

    for ii = 1:row
        cur_seq = Valid_TBPs{iter}(ii,:);        
        Extended_seq = [cur_seq, cur_seq]; 
        for shift = 1:N-1
            cyclic_seq = Extended_seq(1+shift:N+shift);
            key_cyclic_seq = dec2bin(bi2de(cyclic_seq),N);
            if isequal(cyclic_seq, cur_seq) % termination condition for cyclic shift
                break
            end
            if ~isKey(HashTable, key_cyclic_seq)
                Valid_TBPs{iter} = [Valid_TBPs{iter};cyclic_seq]; % find a new TBP
                HashTable(key_cyclic_seq) = 1;
            end
        end
    end
end


TBP_node.list = Valid_TBPs;

% Miscellaneous: Compute total number of TBPs found
aggregate = 0;
for dist = 1:d_tilde
    aggregate = aggregate +size(Valid_TBPs{dist},1);
end

TBP_node.aggregate = aggregate;

% Save results
file_name = ['TBP_node_CC_',code_string,'_N_',num2str(N),'.mat'];

save(file_name,'TBP_node','-v7.3');



