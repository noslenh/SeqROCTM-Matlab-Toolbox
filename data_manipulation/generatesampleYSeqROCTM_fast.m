function [Y, qemp] = generatesampleYSeqROCTM_fast(idx_contexts_X, q, A)
%GENERATESAMPLEYSEQROCTM_FAST Generates the response sequence of a SeqROCTM from a sequence of indexes of the context at each position of the inputs sequence, and the distributions associated to the contexts.
%
% Inputs
%
% idx_contexts  : sequence of indexes of the contexts on each position in
%                 the input sequence
% q             : distributions associated to the contexts
% A             : alphabet
%
% Outputs
%
% Y             : sequence of responses
% qemp          : empirical distributions computed on the simulated samples
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 06/2020

n = length(idx_contexts_X);
Y = -1*ones(1, n);
qemp = zeros(size(q));

next_pos = 1;
% looks for the last context and add the next symbol according to its probability
% distribution
while next_pos < n + 1
    idx_ctx = idx_contexts_X(next_pos);
    [next_symbol, idx] = sampleDiscreteDist(A, q(idx_ctx,:), 1);
    Y(next_pos) = next_symbol;
    next_pos = next_pos + 1;
    qemp(idx_ctx, idx) = qemp(idx_ctx, idx) + 1;        
end
    
% normalize probabilities
qemp = bsxfun(@rdivide,qemp,sum(qemp,2));

end