 function [crc_gen_poly, TBP_node] = Search_DSO_CRC_by_Construction(code_generator, m, d_tilde, N, base, tbp_node)

%
%   The function searches the distance-spectrum-optimal (DSO) CRC generator
%   polynomial for the given TBCC.
%
%   Inputs:
%       1) code_generator: a matrix specifying the generator of TBCC
%       2) m: the degree of the objective DSO CRC generator polynomial
%       3) d_tilde: the distance threshold
%       4) N: a scalar denoting the trellis length
%       5) base: a scalar denoting the final presentation of CRCs,
%       typically 8 or 16.
%       6) tbp_node: the previously cached TBP_node, which can be reused to
%       find other DSO CRCs if valid TBPs remains exactly the same.
%
%   Outputs:
%       1) crc_gen_poly: the DSO CRC generator polynomial in octal
%       2) TBP_node: a struct including the following members:
%           (1) crc_distance: a scalar denoting the minimum undetected
%                   distance provided by the DSO CRC.
%           (1) undetected_spectrum: a 2^(m-1)*d_tilde matrix indicating the 
%                   undetected spectrum for all CRC candidates over distance 
%                   horizon;
%           (2) list: a d_tilde*1 cell indicating the list of length-N TBPs
%                   arranged in increasing distances. The true distance is
%                   the corresponding index minus 1.
%           (3) aggregate: a scalar indicating the number of length-N TBPs
%                   of distance less than 'd_tilde'.
%
%   Notes:
%       1) Must run "Collect_Irreducible_Error_Events" first if IEEEs are
%       not generated before
%       2) The distance index equals true distance + 1, whereas the length
%       index equals true length.
%       

%   Copyright 2020 Hengjie Yang

if nargin <= 4
    base = 8;
end


    

crc_gen_poly = NaN;
TBP_node = {};
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


if nargin < 6 % If no cached TBP_node can be used
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
else
    Valid_TBPs = tbp_node.list;
end
    
%Step 3: Search for DSO CRC generator polynomial
disp('Step 3: Search for the DSO CRC generator polynomial.');

Candidate_CRCs = dec2bin(0:2^(m-1)-1) - '0';
Candidate_CRCs = [ones(2^(m-1),1), Candidate_CRCs, ones(2^(m-1),1)]; % degree order from highest to lowest
Undetected_spectrum = inf(2^(m-1), 1); % each column represents the undected spectrum

Candidate_poly_octal=dec2base(bin2dec(num2str(Candidate_CRCs)),base); % octal form

mask = true(size(Candidate_CRCs,1),1);
locations = find(mask == true);
crc_gen_poly_vec = zeros(1, m+1);


for dist = 2:d_tilde % skip checking all-zero TBPs
    Undetected_spectrum = [Undetected_spectrum, inf(2^(m-1), 1)];
    weight_vec = zeros(size(locations, 1), 1);
    if ~isempty(Valid_TBPs{dist})
        parfor i = 1:size(locations, 1) % This part is parallelizable
            weight_vec(i) = Check_divisible_by_distance(Candidate_CRCs(locations(i),:),Valid_TBPs{dist});
        end

        for i = 1:size(locations,1)
            Undetected_spectrum(locations(i),dist) = weight_vec(i);
        end

        min_weight = min(weight_vec);
        locations = locations(weight_vec == min_weight);
        disp(['    Current distance: ',num2str(dist-1),' number of candidates: ',num2str(size(locations,1))]);
        if length(locations) == 1
            crc_gen_poly = Candidate_poly_octal(locations(1),:);
            crc_gen_poly_vec = Candidate_CRCs(locations(1),:);
            break
        end  
    end
    if dist == d_tilde && length(locations) > 1
        disp(['    d_tilde is insufficient to find the DSO CRC...']);
    end

end



% Step 4: Identify the minimum undetected distance.
disp('Step 4: Identify the minimum undetected distance by the DSO CRC.');
min_dist = -1;
if sum(crc_gen_poly_vec) == 0
    disp('    Failed to find the DSO CRC...');
else
    for dist = 2:d_tilde
        if ~isempty(Valid_TBPs{dist})
            w = Check_divisible_by_distance(crc_gen_poly_vec, Valid_TBPs{dist});
            if w > 0
                min_dist = dist - 1;
                disp(['    DSO CRC polynomial: ',num2str(crc_gen_poly)]);
                disp(['    Minimum undetected distance: ',num2str(min_dist)]);
                break
            end
            if w == 0 && dist == d_tilde
                disp('    d_tilde is insufficient to determine the minimum undetected distance.');
            end
        end
    end
end

TBP_node.crc_distance = min_dist;
TBP_node.undetected_spectrum = Undetected_spectrum;
TBP_node.list = Valid_TBPs;

% Miscellaneous: Compute total number of TBPs found
aggregate = 0;
for dist = 1:d_tilde
    aggregate = aggregate +size(Valid_TBPs{dist},1);
end

TBP_node.aggregate = aggregate;

% Miscellaneous: Compute the state spectrum
v= log2(NumStates);
trellis = poly2trellis(v+1,code_generator);
for dist = 1:d_tilde
    for ii = 1:size(Valid_TBPs{dist},1)
        extended_inputs = [Valid_TBPs{dist}(ii, N-v+1:N), Valid_TBPs{dist}(ii,:)];
        extended_states = zeros(length(extended_inputs),1);
        extended_states(1) = 1;
        for j = 1:length(extended_inputs)-1
            extended_states(j+1) = trellis.nextStates(extended_states(j), extended_inputs(j)+1)+1;
        end
        extended_states = extended_states(v+1:end);
        for i = 1:NumStates
            if any(extended_states == V(i))
                State_spectrum(V(i)) = State_spectrum(V(i)) + 1;
                break
            end
        end
    end
end       


end



function weight = Check_divisible_by_distance(poly_vec,error_events)

% This function computes the undetected weight for "poly_vec" based on "error_events".

weight = 0;
poly_vec = fliplr(poly_vec); % flip degree order from lowest to highest

for i = 1:size(error_events,1)
    [~, remd] = gfdeconv(error_events(i,:),poly_vec, 2);
    if remd == 0
        weight = weight + 1;
    end
end

end

            




        
        
    
    
        
        
        
    












