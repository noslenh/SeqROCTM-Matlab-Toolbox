function count = countctx(contexts, X, alphabet)
%COUNTCTXT gives the number of times and positions where the contexts in [contexts]
%          appears in the sequence X. If the sequence contains a past sequence which 
%          have no context associated to it, it is also returned in count.
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


% Initialize the variable count
ncontexts = length(contexts);
lw = cellfun(@(x) length(x), contexts);
height = max(lw);

T = {};     % past sequences
I = {};     % their positions in the sequence
Tidx = {};  % position in the list of contexts (0 if it is not a context)

if ncontexts ~= 0
    for a = alphabet
        [f, id, c] = find_context(a, alphabet, 2:size(X,2), X, contexts, ncontexts, height);
        T = [T, f];
        I = [I, id];
        Tidx = [Tidx, c];
    end
end

lc = length(Tidx);
count = cell(2, ncontexts + lc); %preallocate for speed
count(1,:) = {0};

% put the results in the structure of count
additional = 0;
for i = 1 : length(Tidx)
    if Tidx{i} > 0
        count{1, Tidx{i}} = length(I{i}); 
        count{2, Tidx{i}} = I{i} + lw(Tidx{i}) - 1;
    else
        position = I{i} + length(T{i}) - 1;
        if position >= height    % add only if the past sequence is large enough to find a context
            additional = additional + 1;
            count{1, ncontexts+additional} = length(I{i});
            count{2, ncontexts+additional} = position;
        end
    end
end
%shrink memory
count(:, ncontexts+additional+1:end) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% version 1 => Aline
% for k = 1 : ncontexts % for each context
%     s = length(contexts{k});
%     for n = 1 : lX-s   % the contexts are count until length(X)-1   
%         if isequal( X(n:n+s-1), contexts{k} ) 
%             count{1,k} = count{1,k} + 1; 
%             count{2,k}(count{1,k}) = n+s-1; 
%         end
%     end
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%% version 2 => Noslen (speed-up)
% 
% lX = length(X);
% count = cell(2, ncontexts);
% if ~isempty(contexts)
%     for n = 1 : ncontexts
%         count{1,n} = 0;             % count ctxt
%         count{2,n} = zeros(1,lX);   % indexes of ctxts (pre-allocating for speed-up)      
%     end
% end
% 
% lc = cellfun(@(x) length(x), contexts);
% max_length = max(lc);
% for i = max_length : lX-1
%     % find the context
%     c = 1;
%     while (c <= ncontexts) && ( ~isequal(contexts{c}, X(i-lc(c)+1:i)) )
%         c = c + 1;
%     end
%     if c <= ncontexts
%         count{1,c} = count{1,c} + 1;
%         count{2,c}(count{1,c}) = i;
%     end
% end
% %shrinking the allocated memory
% for c = 1 : ncontexts
%     count{2,c}(count{1,c}+1:end) = [];
% end    

end