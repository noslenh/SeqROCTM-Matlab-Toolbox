function TEST = getTESTstructure(T, I, la, Y)
%GETTESTSTRUCTURE Gives an structure to be used in the function CTestimator. This function is used only to speed-up
%                 other functions. 
%
% Inputs
%
%   T           : complete tree
%   I           : indexes where the leaves of the complete tree appears in
%                   the sample Y
%   la          : length of the alphabet
%   Y           : sequence of data
%
% Output
%
%   TEST        : structure (contain the complete tree T organized by
%                   branches and levels)
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 07/2020


max_height = max(cellfun(@(x) length(x), T));

% initialize the structure
TEST = cell(max_height+1, 1);
max_level = 0;

% for each element of the complete tree
for i = 1 : length(T)
    % level of that element based on its length
    level = length(T{i}) + 1;
    if level > max_level, max_level = level; end
    % get the statistics corresponding to that leaf (number of times it
    % appears and the transition to the symbols of the alphabet)
    [Nw, Nwa] = get_counts(T{i}, I{i}, Y, la);
    nr = size(TEST{level},2);
    
    % looks for the siblings to put together (in a branch) 
    found = false;
    n = 1;
    while ~found && n <= nr
        if isequal(T{i}(2:end), TEST{level}{1,n}{1,1}(2:end))
            found = true;
            TEST{level}{1,n} = [TEST{level}{1,n}, T{i}];  % add the leave
            TEST{level}{2,n} = [TEST{level}{2,n}, Nw];    % add the counts
            TEST{level}{3,n} = [TEST{level}{3,n}, Nwa];
        else
            n = n + 1;
        end
    end
    if ~found
        TEST{level}{1, nr + 1} = T(i);
        TEST{level}{2, nr + 1} = Nw;
        TEST{level}{3, nr + 1} = Nwa;
    end
end
end

function [Nw, Nwa] = get_counts(w, ind, X, length_alphabet)
    
    Nwa = zeros(length_alphabet,1);
    lw = length(w);
    lx = length(X);
    
    for i = 1 : length(ind)
        if ind(i) + lw <= lx            % this is because ind+l(w) gives the position after w
            loc = X(ind(i) + lw) + 1;   % faster way, interpreting the symbol as index
            Nwa(loc) = Nwa(loc) + 1;
        end  
    end
    
    Nw = sum(Nwa);

end
