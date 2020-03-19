function Poly_node = Search_DSO_CRC(code_generator, d_tilde, N, m, cur_poly_node, base)

%
%   The function serves as a subfunction of
%   "Search_DSO_CRC_by_Construction.m" to search for the
%   distance-specturm-optimal (DSO) CRC polynomial.
%
%   Inputs:
%       1) code_generator: a matrix specifying the generator of TBCC
%       2) d_tilde: the distance threshold
%       3) N: a scalar denoting trellis length
%       4) m: a scalar denoting the CRC degree
%       5) cur_poly_node: existing Poly_node found previously. See format
%               in the Outputs
%       6) base: a scalar denoting the final presentation of CRCs,
%       typically 8 or 16.
%
%   Outputs:
%       1) Poly_node: a struct including the following members:
%           (1) success_flag: True if found the unique DSO CRC, False otherwise.
%           (2) crc_gen_polys: the final list of CRC candidates; the 
%                       list size = 1 (i.e., DSO CRC poly) if success_flag = True, 
%                       and list size >1 otherwise;
%           (3) stopped_distance: -1 if success_flag = True, and a scalar
%                   indicating the stopped distance at which the DSO CRC is
%                   still not finalized if success_flag = False;
%           (4) crc_distance: a positive integer if undetected errors are
%                   found within "d_tilde", and -1 otherise;
%
%   Notes:
%       1) Must run "Reconstruct_TBPs" first if TBPs are not generated before.
%       2) If loaded with existing Poly_node, "cur_poly_node", the algorithm 
%           will extract the stopped distance and refined candidate list and 
%           base on these information to continue the search.
%

%   Copyright 2020 Hengjie Yang

tic

Poly_node = {};
if nargin < 5
   cur_poly_node = {}; 
end

if nargin < 6
    base = 16;
end

code_string = '';
for iter = 1:size(code_generator,2)
    code_string = [code_string, num2str(code_generator(iter)), '_'];
end

file_name = ['TBP_node_TBCC_',code_string,'d_',num2str(d_tilde),'_N_',num2str(N),'.mat'];
if ~exist(file_name, 'file')
    disp(['Error: the file ',file_name, ' does not exist!']);
    return
end


load(file_name, 'TBP_node');

Valid_TBPs = TBP_node.list;
success_flag = false;
stopped_distance = -1;
min_dist = -1;

% Step 3: Search for the DSO CRC generator polynomial.
disp('Step 3: Search for the DSO CRC generator polynomial.');

d_start = 2;
List_size = 2^(m-1);
Candidate_CRCs = dec2bin(0:List_size-1) - '0';
Candidate_CRCs = [ones(List_size,1), Candidate_CRCs, ones(2^(m-1),1)]; % degree order from highest to lowest
Undetected_spectrum = inf(List_size, 1); % each column represents the undected spectrum

Candidate_poly_octal=dec2base(bin2dec(num2str(Candidate_CRCs)),base); % octal form
mask = true(size(Candidate_CRCs,1),1);
locations = find(mask == true);



if ~isempty(cur_poly_node)
    Candidate_poly_octal = cur_poly_node.crc_gen_polys;
    d_start = cur_poly_node.stopped_distance + 1;
    List_size = size(Candidate_poly_octal, 1);
    Candidate_CRCs = dec2bin(base2dec(Candidate_poly_octal, base))-'0';
    Undetected_spectrum = inf(List_size, d_start-1);
    mask = true(size(Candidate_CRCs,1),1);
    locations = find(mask == true);
end




crc_gen_polys = [];
crc_gen_poly_vecs = [];


for dist = d_start:d_tilde % skip checking all-zero TBPs
    Undetected_spectrum = [Undetected_spectrum, inf(List_size, 1)];
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
        disp(['    Current distance: ',num2str(dist-1),' number of candidates: ',...
            num2str(size(locations,1))]);
        if length(locations) == 1
            crc_gen_polys = Candidate_poly_octal(locations(1),:);
            crc_gen_poly_vecs = Candidate_CRCs(locations(1),:);
            success_flag = true;
            break
        end  
    end
    if dist == d_tilde && length(locations) > 1
        crc_gen_polys = Candidate_poly_octal(locations,:);
        crc_gen_poly_vecs = Candidate_CRCs(locations,:);
        stopped_distance = d_tilde;
        disp(['    d_tilde is insufficient to find the DSO CRC... ']);
        disp(['    Stopped distance: ',num2str(stopped_distance),...
            ' # of candidate polynomials: ',num2str(size(crc_gen_polys,1))]);
    end
end



% Step 4: Identify the minimum undetected distance.
if success_flag == true
    disp('Step 4: Identify the minimum undetected distance by the DSO CRC.');
    for dist = 2:d_tilde
        if ~isempty(Valid_TBPs{dist})
            w = Check_divisible_by_distance(crc_gen_poly_vecs(1,:), Valid_TBPs{dist});
            if w > 0
                min_dist = dist - 1;
                disp(['    DSO CRC polynomial: ',num2str(crc_gen_polys(1,:))]);
                disp(['    Minimum undetected distance: ',num2str(min_dist)]);
                break
            end
            if w == 0 && dist == d_tilde
                disp('    d_tilde is insufficient to determine the minimum undetected distance.');
            end
        end
    end
end

% Save results
Poly_node.success_flag = success_flag;
Poly_node.crc_gen_polys = crc_gen_polys;
Poly_node.stopped_distance = stopped_distance;
Poly_node.crc_distance = min_dist;

file_name = ['Poly_node_TBCC_',code_string,'d_',num2str(d_tilde),'_N_',num2str(N),'_m_',num2str(m),'.mat'];
save(file_name,'Poly_node','-v7.3');

timing = toc;
disp(['Execution time: ',num2str(timing),'s']);

end

function weight = Check_divisible_by_distance(poly_vec,error_events)

% This function computes the undetected weight for "poly_vec" based on "error_events".

weight = 0;
poly_vec = fliplr(poly_vec); % flip degree order from lowest to highest

for i = 1:size(error_events,1)
    temp = double(error_events(i,:));
    [~, remd] = gfdeconv(temp,poly_vec, 2);
    if remd == 0
        weight = weight + 1;
    end
end

end



