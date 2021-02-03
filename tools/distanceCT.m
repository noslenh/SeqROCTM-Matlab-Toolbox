function d = distanceCT(tree1, tree2, A)
%DISTANCECT compute a distance between context trees (comparing only the
%           contexts). This distance penalizes the nodes that are in tree1
%           and aren't in tree2 and vice versa according to a weight given
%           to each node.
%
% Inputs   
%
%   tree1 : first context tree 
%   tree2 : second context tree
%       A : alphabet
%
% Outputs
%       d : distance value

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 02/2019

ncontexts1 = length(tree1);
ncontexts2 = length(tree2);

length_alphabet = length(A);

lengths1 = cellfun(@(x) length(x), tree1);
lengths2 = cellfun(@(x) length(x), tree2);
max_height = max([lengths1, lengths2]);

d = 0;

% some auxiliary variables
 powers = length_alphabet.^(0:max_height);
 partial_sum = [0, 0, cumsum(powers(2:end))];
 
 % initialize the weights given to nodes. weights are generated according 
 % to the level of the tree in which a node is
 weights = (length_alphabet-1)./length_alphabet.^(2*(1:max_height) + 1);
 
% vector containing: -1: node exist in both trees
%                     1: node exist only in one of them
%                     0: node does not exist in either of the two trees  
% (for big trees this representation could be sparse)
nodes_xor = zeros(1, partial_sum(end));

% index of the father of each node
idx_father = zeros(1, partial_sum(end));

% initialize nodes_xor with 1 in the nodes that appears in tree1
for c = 1 : ncontexts1
    % initialize node with a context
    node = tree1{c};
    l_node = length(node);
    idx_node = node2index(node, powers(1:l_node), partial_sum(l_node+1));
    % iterate to put the path until the root
    while (~isempty(node)) && (nodes_xor(idx_node) == 0) 
        nodes_xor(idx_node) = 1;
        idx_old = idx_node;
        % sum its weight to the distance (the weight depends only on the height of the tree the node is)
        d = d + weights(l_node);
        % take the father
        node(1) = [];
        l_node = l_node - 1;
        idx_node = node2index(node, powers(1:l_node), partial_sum(l_node+1));
        % store the index of the father
        idx_father(idx_old) = idx_node;
    end   
end

% update the information in nodes_xor using the nodes in tree2: 
% If the node appears in tree1 => nodes_xor = -1
% If the node does not appear in tree1 => nodes_xor = 1
% For the nodes that are not in tree2 
% If it does not appear in tree1 => nodes_xor = 0 (as it was initialized)
% If it appears in tree1 => nodes_xor = 1 (as it was initialized)
for c = 1 : ncontexts2
    % initialize node with a context
    node = tree2{c};
    l_node = length(node);
    idx_node = node2index(node, powers(1:l_node), partial_sum(l_node+1));
 
    finish = false;
    while ~finish && (~isempty(node)) 
        switch nodes_xor(idx_node)
            case 0      % node is not in tree 1         
                nodes_xor(idx_node) = 2;
                % update the distance
                d = d + weights(l_node);
                % take the father to be analyzed (i.e., delete the first symbol)
                node(1) = []; 
                l_node = l_node - 1;
                idx_node = node2index(node, powers(1:l_node), partial_sum(l_node+1));
            case 1      % node is in tree1 
                % put all the path to the root in -1 and update the
                % distance
                root = false;
                while ~root && (nodes_xor(idx_node)~=-1)
                    % mark the node
                    nodes_xor(idx_node) = -1;
                    % update the distance
                    d = d - weights(l_node);
                    % take the index of the father
                    idx_node = idx_father(idx_node);
                    l_node = l_node - 1;
                    if idx_node == 0, root = true; end
                end
                %
                finish = true;
            case 2
                finish = true;
            case -1
                finish = true;
            
        end
    end 
end
end

function idx = node2index(w, powers, ps)

   if isempty(w)
       idx = 0;
   else
       level_idx = sum(w.*powers) + 1;
       idx = level_idx + ps;
   end
end
    