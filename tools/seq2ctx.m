function ctx = seq2ctx(X, tau)
%SEG2CTX Converts a sequence X in the corresponding sequence of indices of contexts
%
% INPUT
%       X   : a sequence of symbols (context tree model)
%       tau : context tree compatible with X
%
% OUTPUT
%
%       ctx : sequence with the indices of the contexts
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 04/2020

max_length = max(cellfun(@(x) length(x), tau));
ctx = -1 * ones(1, length(X)-max_length+1);

for i = max_length : length(X)
    [~, idx] = contextfunction(X(i-max_length+1:i),tau);
    ctx(i-max_length+1) = idx;
end

