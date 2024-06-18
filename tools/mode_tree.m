function leaves = mode_tree(Trees, A)
%MODE_TREE Computes the mode context tree of a set of context tree. 
%
% Inputs
%
%   Trees      : cell array with the context tree used to compute the mode
%   A          : alphabet
%
% Outputs
%
%   leaves     : set of contexts of the mode tree
% 
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 01/2020


    length_alphabet = numel(A);
    
    % maximum height of the trees in Trees
    nTrees = length(Trees);
    
    max_height = 0;
    for t = 1 : nTrees
        tree = Trees{t};
        if ~isempty(tree)
            h = max(cellfun(@(x) length(x), tree, 'uniformOutput', true));
            if  h > max_height, max_height = h; end
        end
    end

    % auxiliary variables
    powers = length_alphabet.^(0:max_height);
    partial_sum = [0, 0, cumsum(powers(2:end))];
    max_num_nodes = sum(powers);

    % variable to store the number of times each node of the complete tree appear as a context
    % (the root have index zero and its frequency is store in the first
    % position of the array)
    freq_contexts = zeros(max_num_nodes, 1);

    for t = 1 : nTrees
        tree = Trees{t};
        ncontexts = length(tree);
        
        if ncontexts == 0   % if the tree is the empty tree, increase the frequency of the root node
            freq_contexts(1) = freq_contexts(1) + 1;
        else
            for c = 1 : ncontexts
                % initialize node with a context
                node = tree{c};
                l_node = length(node);
                
                % get the index of that context in the complete tree
                idx_node = node2index(node, powers(1:l_node), partial_sum(l_node+1));

                % update its frequency as a context
                freq_contexts(idx_node+1) = freq_contexts(idx_node+1) + 1;  
            end
        end
    end
    
    % identify the leaves of the mode tree from the array of frequencies
    leaves = is_leaf_in_mode_tree([], freq_contexts(1), A, freq_contexts, powers, partial_sum, max_height);


end

function T = is_leaf_in_mode_tree(node, freq, Alphabet, frequencies, powers, partial_sum, max_height)
%IS_LEAF_IN_MODE_TREE Recursive function to compute the context of the mode tree from the frequency
%
% Inputs
%
%   node        : sequence of symbols
%   freq        : frequency associated to node
%   alphabet    : alphabet
%   frequencies : frequencies of all nodes in the complete tree
%   power       : auxiliary variable
%   partial_sum : auxiliary variable
%   max_height  : height of the complete tree
%
% Outputs
%
%   T          : set of contexts of the mode tree
% 
    T = {};
    
    if (max_height == length(node))  % if "node" is a leave, and its frequency is not zero, "node" is a context
        if freq ~= 0
            T = node;
        end
    else  
        % if there is no node in the branch induced by "node" with greater
        % or equal frequency than "node", and the frequency of "node" is not zero, then "node" is a context 
        if ~greater_than(node, freq, Alphabet, frequencies, powers, partial_sum, max_height)
            if freq ~= 0
                T = node;
            end
        else
            % ask if the sons of "node" are contexts
            for a = Alphabet
                new_node = [a node];
                l_node = length(new_node);
                idx = node2index(new_node, powers(1:l_node), partial_sum(l_node+1));
                t = is_leaf_in_mode_tree(new_node, frequencies(idx+1), Alphabet, frequencies, powers, partial_sum, max_height); 
                T = [T, t];
            end
        end
    end

end

function found = greater_than(node, freq, Alphabet, frequencies, powers, partial_sum, max_height)
%GREATER_THAN check if there is any node in the branch induced by "node" with greater or equal frequency that "node"
%
% Inputs
%
%   node        : sequence of symbols
%   freq        : frequency associated to node
%   alphabet    : alphabet
%   frequencies : frequencies of all nodes in the complete tree
%   power       : auxiliary variable
%   ps          : auxiliary variable
%
% Outputs
%
%   found       : true if there is a node with greater or equal frequency than
%                 "node", false otherwise
% 
    found = false;    
    l_node = length(node);    

    if l_node < max_height
        found = false;
        a = 1;
        while ~found && a <= numel(Alphabet)
            % check if the sons of node have greater or equal frequency
            new_node = [Alphabet(a) node];
            l_new_node = length(new_node);
            idx_new_node = node2index(new_node, powers(1:l_new_node), partial_sum(l_new_node+1));
            
            if frequencies(idx_new_node + 1) >= freq
                found = true;
            else
                % check in the sons of new_node
                found = greater_than(new_node, freq, Alphabet, frequencies, powers, partial_sum, max_height);
            end
            a = a + 1;
        end
    end
    
end

function idx = node2index(w, powers, ps)
%NODE2INDEX Get the index of the node w in the complete tree when the tree
%           is traversed in level-order/breadth-first search (i.e., we
%           visit every node on a level before going to a lower level)
% Inputs
%
%   w           : sequence of symbols
%   power       : auxiliary variable
%   ps          : auxiliary variable
%
% Outputs
%
%   idx         : global index
%
   if isempty(w)
       idx = 0;
   else
       level_idx = sum(w.*powers) + 1;
       idx = level_idx + ps;
   end
end


   

