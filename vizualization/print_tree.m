function print_tree(tree, A)
%PRINT_TREE Print a context tree in the console
%
% Inputs
%
%   tree : context tree
%   A    : alphabet
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 07/2024

tree_instance = contexts_to_tree(tree, A);
treestr = tree_instance.print();
disp(treestr);

end