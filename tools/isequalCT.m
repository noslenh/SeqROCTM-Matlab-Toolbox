function match = isequalCT(tree1, tree2)
%ISEQUALCT Determine if two context trees are equal.
%          (based only in the comparison of the contexts)
%
% Inputs
%
%   tree1 : first context tree 
%   tree2 : second context tree
%
% Outputs
%
%   match : true if the context tree are equal, false otherwise
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 02/2019

match = false;

l = length(tree1);
if (l == length(tree2))
    
    match = true;
    ctx1 = 1;
    while (match)&&(ctx1 <= l)
        ctx2 = 1;
        found = false;
        while (~found)&&(ctx2 <= l)
            found = isequal(tree1{ctx1}, tree2{ctx2});
            ctx2 = ctx2 + 1;
        end
        match = found;
        ctx1 = ctx1 + 1;
    end 
end

end