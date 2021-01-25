function [tree, P, V, outputs] = bic_WCT(X, Alphabet, max_height, c, df, missing, varargin)
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

fast = false;

switch length(varargin) 
    case 1
        Y = varargin{1};
    case 2
        Y = varargin{1};
        if isempty(Y), Y = X; end
        precomputed_stats = varargin{2};
        fast = true;
    case 0
        Y = X;
end

if fast
    lX_no_nan = precomputed_stats{3};
    penalization_factor = -1 * c * log(lX_no_nan);
    
    % call the estimation algorithm
    [tree, P, V, NODES, STATS] = get_maximizingTree_fast([], length(Alphabet), max_height, penalization_factor, df, precomputed_stats{1}, 0, 0, precomputed_stats{2});
    
    % create the structure 'outputs' with some useful additional outputs
    outputs.nodes = NODES;
    outputs.stats = STATS;
    outputs.nonExistingNodes = precomputed_stats{2};
    outputs.XlengthWithoutNaN = precomputed_stats{3};
else
    lX = length(X);
    
    if missing    %there are missing values
        % get the indexes and total of NaN values
        idx_no_nan = find(~isnan(X));
        
        % initialize the ind_father variable (exclude the positions in which the
        % sequence has NaN values)
        ind_father = idx_no_nan(idx_no_nan > max_height);
        lX_no_nan = length(idx_no_nan);
    else
        lX_no_nan = lX;
        idx_no_nan = 1:lX;
        
        % initialize the ind_father variable
        ind_father = max_height+1:lX;
    end
    
    % initialize the common penalization term
    penalization_factor = -1 * c * log(lX_no_nan);
    
    % estimation algorithm
    [tree, P, ~, V, ~, ~, NODES, STATS, non_existing_nodes] = get_maximizingTree([], length(Alphabet), max_height, ind_father, X, penalization_factor, df, 0, Y);
    
    % create the structure 'outputs' with some useful additional outputs
    outputs.nodes = NODES;
    outputs.stats = STATS;
    outputs.nonExistingNodes = non_existing_nodes;
    outputs.nonanIndexes = idx_no_nan;
    outputs.XlengthWithoutNaN = lX_no_nan;
end