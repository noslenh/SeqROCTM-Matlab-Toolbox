function [Xnew, Y, qemp] = generatesampleYSeqROCTM(X, contexts, q, A)
%GENERATESAMPLEYSEQROCTM Generates the response sequence of a SeqROCTM from
%                        the sequence of inputs, the context tree and the
%                        distrubutions associated to the contexts
%
% Inputs
%
% X        : sequence of inputs
% contexts : set of contexts
% q        : distributions associated to the contexts
% A        : alphabet
%
% Outputs
%
% X        : sequence of inputs
% Y        : sequence of responses
% qemp     : empirical distributions computed on the simulated samples
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 06/2020

    n = length(X);
    % length of the context of maximum length
    max_length = max(cellfun(@(x) length(x), contexts));
    
    % initialize the empirical distributions
    qemp = zeros(size(q));
    % initialize the response sequence
    Y = zeros(1, n - max_length);
    % delete from the input sequence the positions for which it is not
    % possible generate response
    Xnew = X(max_length + 1 : end);
    
    for i = 1 : n - max_length
        % find the context
        [~, idx_ctx] = contextfunction(X(i:i+max_length-1), contexts);
        % use the distribution associated to the context to generate the
        % response
        [Y(i), idx] = samplediscretedist(A, q(idx_ctx,:), 1);
        % update the distriution associated to that context
        qemp(idx_ctx, idx) = qemp(idx_ctx, idx) + 1; 
    end
    
    % normalize probabilities
    qemp = bsxfun(@rdivide,qemp,sum(qemp,2));

end