function [Xnew, Y, qemp] = generatesampleYSeqROCTM(X, contexts, q, A)
%GENERATESAMPLEYSEQROCTM Generates the response sequence of a SeqROCTM from a sequence of inputs, the context tree and the distributions associated to the contexts.
%
% Inputs
%
%   X           : sequence of inputs
%   contexts    : set of contexts
%   q           : distributions indexes by the contexts
%   A           : alphabet
%
% Outputs
%
%   Xnew        : sequence of inputs
%   Y           : sequence of responses
%   qemp        : empirical distributions computed on the simulated samples
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 02/2021

n = length(X);

if isempty(contexts)
    Xnew = X;
    Y = samplediscretedist(A, q, n);
    qemp = histc(Y, A); %counts in same order than the alphabet
else
    % length of the context of maximum length
    max_length = max(cellfun(@(x) length(x), contexts));
    
    % initialize the empirical distributions
    qemp = zeros(size(q));
    
    % initialize the response sequence
    Y = zeros(1, n - max_length);
    
    % delete from the input sequence the positions for which there won´t be
    % a corresponding Y value
    Xnew = X(max_length + 1 : end);
    
    for i = 1 : n - max_length
        % find the context
        [~, idx_ctx] = contextfunction(X(i:i+max_length-1), contexts);
        % use the distribution indexed by the context to generate the
        % response
        if idx_ctx ~= -1
            [Y(i), idx] = samplediscretedist(A, q(idx_ctx,:), 1);
            % update the distribution associated to that context
            qemp(idx_ctx, idx) = qemp(idx_ctx, idx) + 1;
        else
            error('The context tree used to generate the sequence Y is not compatible with the sequence X');
        end
    end    
end

% normalize probabilities
qemp = bsxfun(@rdivide,qemp,sum(qemp,2));

end