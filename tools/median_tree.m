function [leaves, mt_contexts] = median_tree(Trees)
%MEDIAN_TREE Compute the median tree from a set of context trees.
%
% Inputs
%
%   Trees       : list of context trees
%
% Outputs
%
%   leaves      : leaves of the median tree
%   mt_contexts : structure with additional information
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 04/2019

% fix parameters (ToDo: put this as parameters of the function)
max_height = 4;
length_alphabet = 3;
max_num_nodes = 120;     %ToDo: compute this automatically

% some auxiliary variables
powers = length_alphabet.^(0:max_height);
partial_sum = [0, 0, cumsum(powers(2:end))];

% variable to store the number of times each branch appear
sum_nodes = zeros(max_num_nodes, 1);
sum_contexts = zeros(max_num_nodes, 1);

labels_nodes = cell(1, max_num_nodes); % INEFICIENT!!! THINK THE FUNCTION index2node()

nTrees = length(Trees);

for t = 1 : nTrees

 tree = Trees{t};
 ncontexts = length(tree);
 
 nodes_in_the_tree = zeros(max_num_nodes, 1);

 for c = 1 : ncontexts
     
     % initialize node with a context
     node = tree{c};
     l_node = length(node);
     idx_node = node2index(node, powers(1:l_node), partial_sum(l_node+1));
     
     % store the label
     labels_nodes{idx_node} = node;
     
     % update the frequency of the context
     sum_contexts(idx_node) = sum_contexts(idx_node) + 1;
     
     % update the nodes that appears due to that contex
     found = false;
     while ~found %&& ver como poner longitud
         
         % check the node
         nodes_in_the_tree(idx_node) = 1;
         
         % set the node in the father
         node = node(2:end);
         l_node = length(node);
         idx_node = node2index(node, powers(1:l_node), partial_sum(l_node+1));
         
         if (idx_node == 0) || nodes_in_the_tree(idx_node) 
             found = true;
         else
             labels_nodes{idx_node} = node;
         end
     end  
 end
 
 % update the sum for the nodes/branches
 sum_nodes = sum_nodes + nodes_in_the_tree;
 
end

% get the contexts that appears more than half of the time
selected_nodes = sum_nodes > nTrees/2;
mt_nodes = labels_nodes(selected_nodes);

% get from mt_nodes the leaves
ss = cellfun(@(x) length(x), mt_nodes, 'uniformoutput', 1);
max_level = max(ss);

leaves = mt_nodes(ss == max_level);

for l = 1 : max_level-1
    % get the nodes of the current level
    tmp1 = mt_nodes(ss == l); 
    % get the nodes of the next level
    tmp2 = mt_nodes(ss == l+1);
    for n = 1 : length(tmp1)
        %verify if has a son
        found = false;
        j = 1;
        while ~found && j < length(tmp2) 
            if isequal(tmp2{j}(2:end), tmp1{n})
                found = true;
            else
                j = j + 1;
            end
        end
        if ~found, leaves = [leaves, tmp1{n}]; end
    end
end


% get the contexts that appear majority
selected_contexts = find(sum_contexts > 0);
mt_contexts(1,:) = labels_nodes(selected_contexts);
mt_contexts(2,:) = num2cell(sum_contexts(selected_contexts)); 

end
 
 
function idx = node2index(w, powers, ps)

   if isempty(w)
       idx = 0;
   else
       level_idx = sum(w.*powers) + 1;
       idx = level_idx + ps;
   end
end