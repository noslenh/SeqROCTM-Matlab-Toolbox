function print_tree(tree)
%PRINT_TREE Print a context tree in the console
%
% Inputs
%
%   tree : context tree
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 01/2021


treestr = "";

for i = 1 : length(tree)
    % get the context
    nodestr = mat2str(tree{i});
    
    % eliminate '[' , ']' and blanck spaces
    nodestr = strrep(nodestr, '[', '');
    nodestr = strrep(nodestr, ']', '');
    nodestr = strrep(nodestr, ' ', '');
    
    % concatenate
    treestr = strcat(treestr, nodestr);
    treestr = strcat(treestr, " ");
end

if strcmp(treestr,"")
    disp('empty')
else
    disp(treestr);
end

end