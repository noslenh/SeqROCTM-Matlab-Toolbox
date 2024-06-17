function B = bootstrap_blocks(X, renewal_block, seq_length, nB)
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
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 06/2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRICK: Avoid that all bootstrapped samples ended at the same context.
% This could happens when the sequence X is periodic. To avoid that, we will
% generated larger bootstrap samples and, randomly, cut them to the
% specified size.
additional_length = ceil(0.3 * seq_length);
tseq_length = seq_length + additional_length;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


lrb = length(renewal_block);

% allocate memory
idx = zeros(tseq_length,1);

% find the renewal blocks in the sample. If a block appears at position 'i'
% the next block is searched from position 'i + lrp'
lX = length(X);
nrenewals = 1;
i = 1;
while i <= lX - lrb
     if isequal(X(i:i+lrb-1), renewal_block)
        idx(nrenewals) = i;
        nrenewals = nrenewals + 1;
        i = i + lrb;
     else
         i = i + 1;
     end
end

% shrink the allocated memory
idx(nrenewals:end) = [];

% create the blocks using the renewal block
l_idx = length(idx);
blocks = cell(l_idx-1,1);
for i = 1 : l_idx-1
    blocks{i} = X(idx(i):idx(i+1)-1);
end

% create the bootstrap samples
nblocks = length(blocks);
lblocks = cellfun(@(x) length(x), blocks);

%allocate memory
tB = -1*ones(nB, tseq_length + max(lblocks));

for b = 1 : nB
    nseq = 0;
    % draw randomly a block and concatenate until the sequence has the
    % required length
    while nseq < tseq_length
        blk = randi(nblocks);
        tB(b, nseq+1 : nseq+lblocks(blk)) = blocks{blk};
        nseq = nseq + lblocks(blk);
    end
end
%shrink the allocated memory
tB(:,tseq_length+1:end) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRICK - end: Randomly choose where to cut the sequences
B = zeros(nB, seq_length);

rand_positions = randi(additional_length, nB, 1);
for b = 1 : nB
    B(b,:) = tB(b, rand_positions(b)+1 : seq_length+rand_positions(b));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

