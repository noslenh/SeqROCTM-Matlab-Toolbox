function [F, I, log_V, P, Nw, Nwa, NODES, STATS, non_existing_nodes, global_idx] = get_maximizingTree(w, lA, l, ind_father, X, lX, penalization_factor, df, global_idx, Y)
% GET_MAXIMIZINGTREE Recursive function to compute the maximizing tree (see
%                    the definition in Csiszar 2005 IEEE Trans. Inf. Theory)
%
% Inputs
%
%   w                       : node
%   lA                      : length alphabet
%   l                       : maximal height of the tree
%   ind_father              : indexes where the father of the node w appears in the
%                               sequence
%   X                       : sequence
%   lX                      : length of the sequence
%   penalization_factor     : penalization (product of common terms in the penalization term)
%   df                      : type of dregree of freedom (0: |A|-1; 1: number of possible transitions)
%   global_idx              : index of the node visited before node w in the tree of all possible
%                               nodes (including the ones with zero frequency)
%
% Outputs
%
%   F           : contexts
%   log_V       : logarithm of V for the contexts (see definition of V in
%                   the article)
%   NODES       : nodes of the complete tree that were analised
%   STATS       : the values [Phat, ProdV, V, Xi] for each of the analysed
%                   nodes
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 01/2021

    % initialize the variables
    F = {};
    I = {};
    log_V = 0;
    P = [];
    NODES = {};
    STATS = [];
    non_existing_nodes = [];
    
    Nw = 0;
    Nwa = zeros(1,lA);
    
    global_idx = global_idx + 1;
    
    if isempty(w)
%         ind = (l+1 : lX);
        ind = (2 : lX);
    else
        ind = is_in_sample(w, ind_father, X);
    end
    
    if numel(ind) > 0  %~isempty(ind)
              
        if length(w) == l % if w is at maximum level 
            % store w as a context
            F = w;
            I = ind;
            
            % compute Nw and Nwa
            [Nw, Nwa] = get_counts(w, ind, Y, lA);
            P = [P; Nwa/Nw];
            
            % penalization term
            if df == 0
                degree_freedom = (lA - 1); %(lA - 1)/2;
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
            Nw = 0;
            Nwa = zeros(1,lA);
            for a = (0 : lA-1)
                [f, i, logv, p, nw, nwa, nodes, stats, nn, global_idx] = get_maximizingTree([a w], lA, l, ind, X, lX, penalization_factor, df, global_idx, Y);
                F = [F, f];
                I = [I, i];
                P = [P; p];
                non_existing_nodes = [non_existing_nodes; nn];
               
                log_prodV = log_prodV + logv;
                Nw = Nw + nw;
                Nwa = Nwa + nwa;
                % store the statistics
                NODES = [NODES; nodes];
                STATS = [STATS; stats];
            end
            
            % penalization term
            if df == 0
                degree_freedom = (lA - 1); %(lA - 1)/2;
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
                I = ind;
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
    else
        non_existing_nodes = global_idx;
    end
end

function ind = is_in_sample(w, ind_father, X)  % if ind = [], w is not in the sample
    
    lf = length(ind_father);
    ind = zeros(1,lf);
    
    counter = 0;
    for i = 1 : lf
        idxson = ind_father(i) - 1;
        if (idxson > 0) && (w(1) == X(idxson))
            counter = counter + 1;
            % index at which w appears in X (the beginning of w)
            ind(counter) = idxson;
        end
    end
    % shrink the allocated memory
    ind(counter+1:end) = [];
end

function [Nw, Nwa] = get_counts(w, ind, X, length_alphabet)
%GET_COUNTS Gives how often w appears in X and how often each symbol of the
%           alphabet appears after w.
%
% Inputs
%   w               : sub-sequence
%   ind             : position where w happens in the sequence
%   X               : sequence
%   length_alphabet : length of the alphabet
%
% Outputs
%   Nw              : Number of occurrences of w in X
%   Nwa             : Number of occurrences of each symbol in the alphabet
%                     after w 

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 05/2019

    Nwa = zeros(1, length_alphabet);
    lw = length(w);
    lx = length(X);
    
    for i = 1 : length(ind)
        pos = ind(i) + lw;
        if pos <= lx            % this is because ind+l(w) gives the position after w
            loc = X(pos) + 1;   % faster way: interpreting the symbol as index
            Nwa(loc) = Nwa(loc) + 1;
        end  
    end   
    Nw = sum(Nwa);
end
