function d = distance_balding(tree1, tree2, all_nodes) 
%DISTANCE Compute the distance of Balding between two context trees. 
%
% Inputs
%
%	tree1 		: first context tree
%	tree2 		: second context tree
%	all_nodes 	: nodes of the complete tree
%	all_weights : weights associated to the nodes of the complete tree
%
% Outputs
%
%	d			: value of the distance
%
% Reference:
%   Balding et. al. 
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr)
%Date   : 10/2020

% compute the weights corresponding to the nodes of the complete tree
n_nodes = length(all_nodes);
level_nodes = cellfun(@length, all_nodes);
max_level = max(level_nodes);

amount_by_level = zeros(1, max_level);
for n = 1 : n_nodes
    amount_by_level(level_nodes(n)) = amount_by_level(level_nodes(n)) + 1;
end

normalization = sum(1 ./ 2.^(1:max_level));
all_weights = (1./(2.^level_nodes) .* 1./amount_by_level(level_nodes)) / normalization;

% initialize the distance value
d = 0;

% context tree to nodes
nodes1 = tree2nodes(tree1);
nodes2 = tree2nodes(tree2);

% number of nodes
l_nodes1 = length(nodes1);
l_nodes2 = length(nodes2);

% variable to store the nodes that are in tree1 but not in tree2, and the nodes that are in tree2 and in tree1, respectively
nodes_in_1_not_in_2 = [];
nodes_in_2_in_1 = [];

% find the nodes that are in tree1 and in tree2 and the nodes that are in tree1 and not in tree2
for n = 1 : l_nodes1
    found = false;
    count = 0;
    while ~found && count < l_nodes2
        count = count + 1;
        if isequal(nodes1{n}, nodes2{count})
            found = true;
            nodes_in_2_in_1 = [nodes_in_2_in_1; count];
        end
    end
    if ~found 
        nodes_in_1_not_in_2 = [nodes_in_1_not_in_2; n];
    end    
end

% find the nodes that are in tree2 and not in tree1
nodes_in_2_not_in_1 = setdiff(1:l_nodes2, nodes_in_2_in_1);

% add the weights corresponding to the nodes that are in tree1 and not in tree2 
for j = 1 : length(nodes_in_1_not_in_2)
    idx = nodes_in_1_not_in_2(j);
    d = d + get_weight(nodes1{idx}, all_nodes, all_weights);
end

% add the weights corresponding 
for j = 1 : length(nodes_in_2_not_in_1)
    idx = nodes_in_2_not_in_1(j);
    d = d + get_weight(nodes2{idx}, all_nodes, all_weights);
end
end

function w = get_weight(node, all_nodes, all_weights)

    found = false;
    count = 0;
    while ~found
        count = count + 1;
        if isequal(node, all_nodes{count})
            found = true;
        end
    end
    w = all_weights(count);
end