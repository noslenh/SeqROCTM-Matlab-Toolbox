function nodes = tree2nodes(tree)
%TREE2NODES Gives the set of nodes corresponding to a context tree
%
% Inputs
%   tree   : a context tree
%
% Outputs
%   nodes  : the set of nodes
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 02/2021

nodes = {};

%number of contexts
nc = length(tree);

for i = 1 : nc
    %get the context, add it to the set of nodes and add the nodes from it
    %to the root
    node = tree{i};
    while ~isempty(node)
        nodes = [nodes, node];
        node = node(2:end);
    end
end

%delete the repeated nodes
N = cellfun(@(x) num2str(x(:)'), nodes,'UniformOutput', false);
[~, idx] = unique(N);
nodes = nodes(idx);

end
