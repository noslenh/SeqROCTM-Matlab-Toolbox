function [tree, P, V, results] = bic_WCT(X, Alphabet, max_height, c, df, missing, varargin)
%BIC_WCT Estimate a context tree model from a sequence X using the BIC
%        criterion
%
% Inputs
%   X           : sequence of symbols taking values in the alphabet
%   Alphabet    : alphabet
%   height      : height of the complete tree
%   c           : penalization constant of the BIC criteria
%   df          : type of degree of freedom function
%   missing     : 1 if treatment of missing values is needed, 0 otherwise
%   varargin    : {1}-> Y, {2}-> precomputed_stats
%
% Outputs
%   tree        : context tree estimated
%   P           : distributions associated to the contexts
%   V           : log(V) values for the contexts (see [2])
%   results     : structure with the following fields:
%       'nodes'             --  all possible nodes of the tree 
%       'stats'             --  statistics of the nodes [log_Ps, \sum log_Vas, log_Vs, \Xi_s, Nw, Nwa]
%       'nonExistingNodes'  --  nodes of the tree that do not appear in the
%                               sequence X
%       'XlengthWithoutNaN' --  length of the sequence X discarding the Nan
%                               values
%       'nonanIndexes'      --  indexes of the non Nan elements in the
%                               sequence X

%   References:  
%      [1] I. Csiszar et al., IEEE Trans. Inform. Theory, 3, 52, 1007-1016 (2006)
%      [2] N. Hernández et al., arXiv:2009.06371, (2021). 

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 01/2021

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
    results.nodes = NODES;
    results.stats = STATS;
    results.nonExistingNodes = precomputed_stats{2};
    results.XlengthWithoutNaN = precomputed_stats{3};
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
    results.nodes = NODES;
    results.stats = STATS;
    results.nonExistingNodes = non_existing_nodes;
    results.nonanIndexes = idx_no_nan;
    results.XlengthWithoutNaN = lX_no_nan;
end