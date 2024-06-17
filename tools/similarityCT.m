function s = similarityCT(tree_ref, tree, weights_ref)
%SIMILARITYCT Compute a similarity between the two context trees
%             tree_reference and tree (comparing only the contexts) by
%             counting which context of tree_reference appears in tree and
%             penalizing by a corresponding weight. 
%
% Inputs   
%
%   tree_ref    : reference context tree 
%   tree        : second context tree
%   weights_ref : weights for the contexts in tree_ref
%
% Outputs
%             s : similarity value
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 10/2020

s = 0;
nc = length(tree_ref);
nc_tree =  length(tree);
for i = 1 : nc
    % check if the ith context appears in tree
    found = false;
    count = 1;
    while ~found && count <= nc_tree
        if isequal(tree_ref{i}, tree{count})
            found = true;
            s = s + weights_ref(i);
        else
            count = count+1;
        end
    end
end

end

