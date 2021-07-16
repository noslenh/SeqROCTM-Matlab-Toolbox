function [BX, BY] = bootstrap_blocks2(X, Y, renewal_block, seq_length, nB)
%BOOTSTRAP_BLOCKS2 Create nB bootstrap samples from the seqROCTM (X,Y)
%                  using a renewal block/context of X. 
%
% Inputs
%
%   X             : input sequence  
%   Y             : response sequence
%   renewal_block : renewal block of sequence X
%   seq_length    : length of the bootstrap samples
%   nB            : number of bootstrap samples
%
% Outputs
%
%   BX            : matrix containing on each row a bootstrap sample of X
%   BY            : matrix containing on each row a corresponding bootstrap
%                   sample of Y 
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 06/2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRICK: Avoid that all bootstrapped samples X ended at the same context.
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
% the next block is searched from position 'i + lrb'
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

% create the blocks using the renewal point
l_idx = length(idx);
blocks = cell(l_idx-1,2);
for i = 1 : l_idx-1
    blocks{i,1} = X(idx(i):idx(i+1)-1);
    blocks{i,2} = Y(idx(i)+lrb : idx(i+1)+lrb-1);
end

% create the bootstrap samples
nblocks = length(blocks);
lblocks = cellfun(@(x) length(x), blocks(:,1));

%allocate memory
tBX = -1*ones(nB, tseq_length + lrb + max(lblocks));
tBY = -1*ones(nB, tseq_length + lrb + 2*max(lblocks));

for b = 1 : nB
    nseq = 0;
    % draw randomly a block and concatenate until the sequence has the
    % required length
    while nseq < tseq_length + lrb
        blk = randi(nblocks);
        tBX(b, nseq+1 : nseq+lblocks(blk)) = blocks{blk,1};
        tBY(b, nseq+1+lrb : nseq+lblocks(blk)+lrb) = blocks{blk,2};
        nseq = nseq + lblocks(blk);
    end
end
%shrink the allocated memory. Also delete the first lrp values because
%those values were not generated for the Y sequences.
tBX(:,1:lrb) = [];
tBX(:,tseq_length+1:end) = [];

tBY(:,1:lrb) = [];
tBY(:,tseq_length+1:end) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRICK - end: Randomly choose where to cut the sequences
BX = zeros(nB, seq_length);
BY = zeros(nB, seq_length);

rand_positions = randi(additional_length, nB, 1);
for b = 1 : nB
    BX(b,:) = tBX(b, rand_positions(b)+1 : seq_length+rand_positions(b));
    BY(b,:) = tBY(b, rand_positions(b)+1 : seq_length+rand_positions(b));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

