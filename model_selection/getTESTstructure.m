function TEST = getTESTstructure(T, I, max_height, la, Y)

    %la = length(Alphabet);
    
    TEST = cell(max_height+1,1);
    max_level = 0;
    
     for i = 1 : length(T)
        
        level = length(T{i}) + 1;
        if level > max_level, max_level = level; end
        
        [Nw, Nwa] = get_counts(T{i}, I{i}, Y, la);
        nr = size(TEST{level},2);
        
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
        if ind(i) + lw <= lx    % this is because ind+l(w) gives the position after w
            % loc = X(ind(i) + length(w)) == Alphabet;
            loc = X(ind(i) + lw) + 1; % faster way, interpreting the symbol as index
            Nwa(loc) = Nwa(loc) + 1;
        end  
    end
    
    Nw = sum(Nwa);

end
