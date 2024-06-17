function [F, P, log_V, NODES, STATS, global_idx, local_idx] = get_maximizingTree_fast(w, lA, max_height, penalization_factor, df, precomputed_stats, ...
                                                                            global_idx, local_idx, non_existing_nodes)
% GET_MAXIMIZINGTREE_FAST A faster version of the recursive function GET_MAXIMIZINGTREE
%                         in which the statistics associated to the nodes are given.  
%
% Inputs
%
%   w                       : node
%   lA                      : length alphabet
%   l                       : maximal height of the tree
%   lX                      : length of the sequence
%   penalization_factor     : penalization (product of common terms in the penalization term)
%   df                      : type of degree of freedom (0: |A|-1; 1: number of possible transitions)
%   precomputed_stats       : statistics of all the nodes of the complete
%                               tree to speed-up computations
%   global_idx              : index of the node visited before node w in the tree of all possible
%                               nodes (including the ones with zero frequency)
%   local_idx               : index of the node visited before node w in the complete tree
%   non_existing_nodes      : global index of the nodes with zero frequency
%
% Outputs
%
%   F                       : contexts
%   log_V                   : logarithm of V for the contexts (see definition of V in
%                               the article)
%   NODES                   : nodes of the complete tree that were analyzed
%   STATS                   : the values [Phat, ProdV, V, Xi] for each of the analyzed
%                               nodes
%   global_idx              : index of w in the tree of all possible nodes (including
%                               the ones with zero frequency)
%   local_idx               : index of w in the complete tree
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 01/2021

    % initialize the variables
    F = {};
    log_V = 0;
    P = [];
    NODES = {};
    STATS = [];
    
    % update the global marker of the node (this marker count even
    % zero-frequency nodes)
    global_idx = global_idx + 1;
        
    if ~ismember(global_idx, non_existing_nodes)
        
        if length(w) == max_height % if w is at maximum level
            
            % update the local marker (this marker does not count zero
            % frequency nodes). It is used to known the position of the
            % node in the precomputed_stats matrix
            local_idx = local_idx + 1;
            
            % store w as a context
            F = w;
            
            % get Nw and Nwa from the matrix
            Nw = precomputed_stats(local_idx, 1);
            Nwa = precomputed_stats(local_idx, 2 : lA+1);
            
            % compute the probs
            P = [P; Nwa/Nw];
            
            % penalization term
            if df == 0
                degree_freedom = (lA - 1);
            else
                degree_freedom = sum(Nwa > 0);
            end
            log_flag = penalization_factor * degree_freedom;
            
            % compute log(V)
            idxp = Nwa > 0;
            log_V = log_flag + sum( Nwa(idxp).* log(Nwa(idxp)/Nw) );
            
            % store the statistics
            NODES = w;
            STATS(1:3) = log_V;
            STATS(4) = 0;
            STATS(5) = Nw;
            STATS(6:6+lA-1) = Nwa;
        else
            
            % compute the sum of log(L) for the sons of w
            log_prodV = 0;
  
            for a = (0 : lA-1)
                [f, p, logv, nodes, stats, global_idx, local_idx] = get_maximizingTree_fast([a w], lA, max_height, penalization_factor, df, ...
                                                                            precomputed_stats, global_idx, local_idx, non_existing_nodes);
                F = [F, f];
                P = [P; p];
                
                log_prodV = log_prodV + logv;
                
                % store the statistics
                NODES = [NODES; nodes];
                STATS = [STATS; stats];
            end
            
            local_idx = local_idx + 1;
            
            % get Nw and Nwa from the matrix
            Nw = precomputed_stats(local_idx, 1);
            Nwa = precomputed_stats(local_idx, 2 : lA+1);
            
            % penalization term
            if df == 0
                degree_freedom = (lA - 1);
            else
                degree_freedom = sum(Nwa > 0);
            end
            log_flag = penalization_factor * degree_freedom;
            
            % compute log(L) for the tree with w
            idxp = Nwa > 0;
            log_L = log_flag + sum( Nwa(idxp).* log(Nwa(idxp)/Nw) );
            
            if isempty(F)||(log_prodV <= log_L) % X = 0, discard previous contexts, new context is w
                % None of the children is leaf, so w is leaf
                log_V = log_L;
                F = w;
                P = Nwa/Nw;
                
                % store the statistics
                STATS = [STATS; log_L log_prodV log_V 0 Nw Nwa];
            else
                % X = 1, keep the previous contexts
                log_V = log_prodV;
                
                % store the statistics
                STATS = [STATS; log_L log_prodV log_V 1 Nw Nwa];
            end
            % store the statistics
            NODES = [NODES; w];
        end
    end