function all_trees = get_search_space_trees(complete_tree, A)
%GET_SEARCH_SPACE_TREES compute all the context trees that could result
%                       from the estimation process (search space). 
% Inputs
%
%   complete_tree : all nodes of the complete tree
%   A                   : alphabet
%
% Output
%
%  all_trees : cell array with all possible context trees
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 02/2021

 % get the nodes of the complete tree
 nodes_complete_tree = tree2nodes(complete_tree);
 
 % compute the height of the complete tree
 max_height = max(cellfun(@length, nodes_complete_tree));   
 
 % call the recursive function at the empty tree
 all_trees = get_trees([], max_height, nodes_complete_tree, A);

end

function set_trees = get_trees(tree, max_height, nodes_complete_tree, A)
%GET_TREES Get all possible context tree that can be constructed by
%          branching the leaves of a tree at the maximum level
%
% Input
%
%  tree                : context tree 
%  max_height          : maximum height to grow 
%  nodes_complete_tree : nodes of the complete tree
%  A                   : alphabet
%
% Output
%
%  set_trees           : cell array with all the context trees that can be
%                        obtained

    % initialize
    set_trees = cell(0);
    
    if isempty(tree)  % if tree is the empty tree (root)
        % branch it to the next level
        new_tree = branch(tree, A, nodes_complete_tree);
        % store the empty tree and the subtree branched from it
        set_trees = {tree; new_tree};
        % call the branch function in the new tree
        set_trees = [set_trees; get_trees(new_tree, max_height, nodes_complete_tree, A)];
    else
        
        % get the height of the current tree
        l_leaves = cellfun(@length, tree);
        heigth_tree = max(l_leaves);
        
        % if the height if the maximum height stop branching
        if heigth_tree == max_height
            set_trees = {};
        else
            
            % looks for the leaves to grow (the leaves at the highest level)
            idx_leaves_lk = find(l_leaves == max(l_leaves));

            % branch each of the leaves
            n_lk = length(idx_leaves_lk);
            sub_trees = cell(n_lk,1);
            %indexes of the leaves that will grow
            idx_possible_branch = zeros(n_lk,1);
            n_new_subtrees = 0;
            for i = 1 : n_lk
                tmp = branch(tree{idx_leaves_lk(i)}, A, nodes_complete_tree);
                % store the new subtree only if it is not empty 
                if ~isempty(tmp)
                    n_new_subtrees = n_new_subtrees + 1;
                    sub_trees{n_new_subtrees} = tmp;
                    idx_possible_branch(n_new_subtrees) = idx_leaves_lk(i);
                end
            end

            % build trees by growing all the combination of the leaves that
            % can grow
            for i = 1 : n_new_subtrees
                % all possible combinations of i subtrees
                leaves_to_branch = nchoosek(1:n_new_subtrees,i);
                for r = 1 : size(leaves_to_branch,1)
                    new_tree = tree;
                    % delete the leaves
                    new_tree(idx_possible_branch(leaves_to_branch(r,:))) = [];
                    % add the subtrees growing from those leaves
                    tmp = sub_trees(leaves_to_branch(r,:));
                    new_tree = [new_tree, tmp{:}];
                    set_trees{end+1,1} = new_tree;
                    % get the trees that can be generated from the new tree
                    tmp = get_trees(new_tree, max_height, nodes_complete_tree, A);
                    % update the final cell array of trees
                    set_trees = vertcat(set_trees, tmp);
                end
            end
        end
    end

end

function new_leaves = branch(leaf, A, nodes_complete_tree)
%BRANCH Return the leaves that can be branched from a given leaf
%
% Inputs
%   leaf                : leaf from which the new leaves will branch
%   A                   : alphabet
%   nodes_complete_tree : nodes of the complete tree
%
% Output
%
%   new_leaves          : cell array with the resulting leaves
%

    new_leaves = {};
    for a = A
        new_leaf = [a leaf];
        % check if the new leaf is possible in the complete tree
        if exist_node(new_leaf, nodes_complete_tree)
            new_leaves = [new_leaves, new_leaf];
        end
    end
end

function exist = exist_node(new_leaf, nodes_complete_tree)
%EXIST_NODE Check if a node exist in the complete tree
%
% Inputs
%   new_leaf            : node to be checked 
%   nodes_complete_tree : nodes of the complete tree
%
% Output
%
%   exist               : true if new_leaf is among the nodes of the
%                           complete tree

 exist = false;
 i = 1;
 while ~exist && i <= length(nodes_complete_tree)
     if isequal(new_leaf, nodes_complete_tree{i})
         exist = true;
     else
         i = i + 1;
     end
 end

end
