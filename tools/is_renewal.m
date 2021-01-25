function rp = is_renewal(pastlength, idw, T, height, nA, lcontexts)
%IS_RENEWAL Recursive function to check if a context is a renewal point

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 01/2021

rp = true;

%if pastlength is greater or equal than height it is always possible to
%generate next symbol
if pastlength < height
    s = 1;
    while rp && s <= nA
        
        %take the possible new context(s)
        idw_new = T{idw, s};
        nc = numel(idw_new);
        
        if nc > 1  %if the new context is not a suffix of idw + symbol
            
            %check if the length of the past is enough to recover one of
            %the possible contexts
            j = 1;
            while rp && j <= nc
                rp = lcontexts(idw_new(j)) <= pastlength;
                j = j + 1;
            end
            %continue checking (if rp is still true)
            j = 1;
            while rp && j <= nc
                rp = is_renewal(pastlength+1, idw_new(j), T, height, nA, lcontexts);
                j = j + 1;
            end
        else %if the new context is a suffix, check them
             rp = is_renewal(pastlength+1, idw_new, T, height, nA, lcontexts);
        end
        s = s + 1;
    end
end

end

