function B = bootstrap_blocks(X, renewal_point, seq_length, nB)
%BOOTSTRAP_BLOCKS Create nB bootstrap samples from a sequence X and a renewal point 
%
% Inputs
%
%   X             : sequence  
%   renewal_point : renewal point
%   seq_length    : length of the bootstrap samples
%   nB            : number of bootstrap samples
%
% Outputs
%
%   B             : matrix containing on each row a bootstrap sample
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 05/2020


lrp = length(renewal_point);

% allocate memory
idx = zeros(seq_length,1);
% find the renewal points in the sample
nrenewals = 1;
for i = 1 : seq_length - lrp
    if isequal(X(i:i+lrp-1), renewal_point)
        idx(nrenewals) = i;
        nrenewals = nrenewals + 1;
    end
end
% shrink the allocated memory
idx(nrenewals:end) = [];

% create the blocks using the renewal point
l_idx = length(idx);
blocks = cell(l_idx-1,1);
for i = 1 : l_idx-1
    blocks{i} = X(idx(i):idx(i+1)-1);
end

% create the bootstrap samples
nblocks = length(blocks);
lblocks = cellfun(@(x) length(x), blocks);

%allocate memory
B = -1*ones(nB, seq_length + max(lblocks));

for b = 1 : nB
    nseq = 0;
    % draw randomly a block and concatenate until the sequence has the
    % required length
    while nseq < seq_length
        blk = randi(nblocks);
        B(b, nseq+1 : nseq+lblocks(blk)) = blocks{blk};
        nseq = nseq + lblocks(blk);
    end
end
%shrink the allocated memory
B(:,seq_length+1:end) = [];
