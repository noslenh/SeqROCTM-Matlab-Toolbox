function [F, I, Fidx] = find_context(w, Alphabet, ind_father, X, contexts, ncontexts, max_height)
%FIND_CONTEXT recursive function to 
% Inputs
%
%   w          : sequence of symbols
%   alphabet   : alphabet
%   max_height : height of the complete tree
%   ind_father : indexes where the father of w appears in the sequence X
%   X          : sequence of data
%
% Outputs
%
%   F          : set of contexts of the complete tree
%   I          : indexes indicating the position of the contexts of the complete
%                  tree in the sequence X
%   Iidx       : total number of pairs of siblings in the complete tree (useful
%                  when the prune is based on statistical testing)
% 
% Usage
%
% 

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 05/2020

    F = {};
    I = {};
    Fidx = {};
    
    ind = is_in_sample(w, ind_father, X);
    if numel(ind) > 0
        [d, c] = is_context(w, contexts, ncontexts);
        if d||(length(w) == max_height)  % if it is a context
            F = w;
            I = ind;
            Fidx = c;
        else
            for a = Alphabet
                [f, i, fidx] = find_context([a w], Alphabet, ind, X, contexts, ncontexts, max_height);
                F = [F, f];
                I = [I, i];
                Fidx = [Fidx, fidx];
            end
            if isempty(F)  % if non of my soon appears, w is a leaf
                F = w;
                I = ind;
                Fidx = c;
            end
        end 
    end
end

function ind = is_in_sample(w, ind_father, X)  % if ind = [], w is not in the sample

    % allocate memory for speed
    lf = length(ind_father);
    ind = zeros(1,lf);
    
    %
    counter = 0;
    for i = 1 : lf
        ii = ind_father(i) - 1;
        if (ii > 0) && (w(1) == X(ii))
            counter = counter + 1;
            ind(counter) = ii;
        end
    end
    % shrink the allocated memory
    ind(counter+1:end) = [];

end

function [d, c] = is_context(w, contexts, ncontexts)
%IS_CONTEXT check if w belong to the set of contexts

% d : true if w belong to contexts
% c : index of the context in the list (0 if d=false)

    d = false;
    c = 1;
    while (c <= ncontexts) && ( ~isequal(contexts{c}, w) )
        c = c + 1;
    end
    if c <= ncontexts
        d = true;
    else
        c = 0;
    end
end