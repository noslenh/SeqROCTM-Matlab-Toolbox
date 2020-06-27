function [Y, sP] = generateSeqYfromContextTree(idx_contexts_X, P, A)
% This function generates a sequence Y according to a realization X of a 
% context tree tau0. For computational reasons, instead of given X and tau0 as input, the
% ordered list of context´s indexes in X is given. 

% the results are the sequence and the sampled transition probabilities matrix

lengthSeq = length(idx_contexts_X);
Y = -1*ones(1, lengthSeq);
sP = zeros(size(P));

if isempty(idx_contexts_X)      % generate an iid sequence (empty tree)
    Y = sampleDiscreteDist(A, P, lengthSeq);
    sP = hist(Y, unique(Y));    %counts in same order than the alphabet
else
    next_pos = 1;
    % looks for the last context and add the next symbol according to its probs
    % distribution
    while next_pos < lengthSeq + 1
        idx_ctx = idx_contexts_X(next_pos);
        [next_symbol, idx] = sampleDiscreteDist(A, P(idx_ctx,:), 1);
        Y(next_pos) = next_symbol;
        next_pos = next_pos + 1;
        sP(idx_ctx, idx) = sP(idx_ctx, idx) + 1;        
    end
    
end
% normalize probabilities
sP = bsxfun(@rdivide,sP,sum(sP,2));
end