function T = admissiblecontexttree(X, max_height, alphabet)
%COMPLETETREE  Compute the complete tree of height max_height compatible
%               with the data X
% Inputs
%
%   X           : sequence of symbols taking values in the alphabet
%   max_height  : height of the complete tree
%   alphabet    : alphabet 
%
% Outputs
%
%   T           : admissible context tree
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 05/2021

% get the complete tree
T = completetree(X, max_height, alphabet);

% delete the contexts that have no sibling
ncontexts = length(T);

for i = 1 : ncontexts
    current_context = T{i};
    % check if current_context has sibling
    found = false;
    j = 1;
    suffix = current_context(2:end);
    lsuffix = length(suffix);
    while ~found 
        if (j~=i) && (length(T{j}) >= lsuffix) && isequal(suffix, T{j}(end-length(suffix)+1:end))
            found = true;
        else
            if j < ncontexts
                j = j + 1; % continue the comparison
            else
                % no sibling was found
                T{i} = T{i}(2:end); %prune
                % initialize all the variables to check the new context
                current_context = T{i};
                j = 1;
                suffix = current_context(2:end);
                lsuffix = length(suffix);
            end
        end
    end
end

