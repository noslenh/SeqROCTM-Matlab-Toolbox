function count = countstr(contexts, X)
% COUNTSTR counts the number of times and positions where the strings in [contexts]
%          appears in the sequence X. 
%
% Inputs
%
%   X        : sequence
%   contexts : set of contexts
%   alphabet : alphabet
%
% Output
% 
%   count    : cell array containing in the first row the number of times
%              the contexts appear in the sequence. In the second row the
%              positions and in the 

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 05/2020


ncontexts = length(contexts);
lX = length(X);

count = cell(2, ncontexts); 	%preallocate for speed
% count(1,:) = {0};
for n = 1 : ncontexts
    count{1,n} = 0;             % count ctxt
    count{2,n} = zeros(1,lX);   % indices of ctxts (preallocating for speed-up)      
end

%%%%%%%%%% version 1 => Aline %%%%%%%%%%
for k = 1 : ncontexts % for each context
    s = length(contexts{k});
    for n = 1 : lX-s   % the contexts are count until length(X)-1   
        if isequal( X(n:n+s-1), contexts{k} ) 
            count{1,k} = count{1,k} + 1; 
            count{2,k}(count{1,k}) = n+s-1; 
        end
    end
end

%shrinking the allocated memory
for c = 1 : ncontexts
    count{2,c}(count{1,c}+1:end) = [];
end    

end