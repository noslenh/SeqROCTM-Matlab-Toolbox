function [F, I, nT] = is_leaf(w, Alphabet, max_height, ind_father, X)
%IS_LEAF recursive function to compute the complete tree
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
%   nT         : total number of pairs of siblings in the complete tree (useful
%                  when the prune is based on statistical testing)
% 
% Usage
%
% 

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 05/2020

    F = {};
    I = {};
	nT = 0;
	nson = 0;
    
    ind = is_in_sample(w, ind_father, X);
    if numel(ind) > 0      % (do not change here the threshold to filter the context by its number of occurrences!!)
        if length(w) == max_height  % if the level is max_height, w is a leaf
            F = w;
            I = ind;
			nT = 0;
        else
            for a = Alphabet
                [f, i, nt] = is_leaf([a w], Alphabet, max_height, ind, X);
                F = [F, f];
                I = [I, i];
				% counting the sons of w
				if ~isempty(f), nson = nson + 1; end
				nT = nT + nt;
            end
            if isempty(F)  % if non of my soon appears, w is a leaf
                F = w;
                I = ind;
				nT = 0;
            end
        end 
        if nson > 1
            nT = nT + nchoosek(nson,2); % update the number of pairs of siblings (test) given the number of sons
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