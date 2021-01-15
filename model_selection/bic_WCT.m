function [tree, P, V, NODES, STATS] = bic_WCT(X, Alphabet, max_height, c, df, varargin)
%BIC_WCT Estimate a context tree model from a sequence using the BIC
%        criterion introduced in Csiszar 2005 IEEE Trans. Inf. Theory 
%
% Inputs
%   X           : sequence of symbols taking values in the alphabet
%   Alphabet    : alphabet
%   height      : height of the complete tree
%   c           : penalization constant of the BIC criteria
%   Y, precomputed_stats
%
% Outputs
%   tree        : context tree estimated
%   idx         : indexes of the context in the sample X
%   V           : log(V) values for the contexts (see the article)
%   NODES       : nodes of the complete tree that were analised
%   STATS       : the values [Phat, ProdV, V, Xi] for each of the analysed
%                   nodes
%  

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 10/2020

fast = true;

switch length(varargin) 
    case 1
        Y = varargin{1};
        fast = false;
    case 2
        Y = varargin{1};
        if isempty(Y), Y = X; end
        precomputed_stats = varargin{2};
    case 0
        Y = X;
        fast = false;
end

lX = length(X);
penalization_factor = -1 * c * log(lX);

if fast
    [tree, V, P, NODES, STATS] = get_maximizingTree_fast([], length(Alphabet), max_height, lX, penalization_factor, df, precomputed_stats{1}, 0, 0, precomputed_stats{2});
else
    [tree, ~, V, P, ~, ~, NODES, STATS] = get_maximizingTree([], length(Alphabet), max_height, [], X, lX, penalization_factor, df, 0, Y);
end



