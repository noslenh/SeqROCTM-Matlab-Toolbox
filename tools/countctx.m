function [count, Nwa] = countctx(contexts, X, A, idx_without_NaN, Y)
%COUNTCTXT Gives the number of times and positions where the contexts appear in a sequence.
%          If the sequence X contains a past sequence which have no context
%          associated to it, it is also returned in count. Besides the
%          number of transitions from each context to each symbol of the
%          alphabet is returned.
%
% Inputs
%
%   contexts            : context tree
%   X                   : sequence from which the contexts are goint to be
%                         counted
%   A                   : alphabet
%   idx_without_NaN     : positions of no Nan values in X (if no value is
%                         given this variable is initialize with all
%                         possible positions)
%   Y                   : sequence from which the transitions are going to
%                         be counted.
%
% Output
% 
%   count    : cell array containing in the first row the number of times
%              each context appear in the sequence. In the second row the
%              positions in X where the context appears. 
%   Nwa      : a matrix containing in each row the transitions from a
%              context to each symbol of the alphabet. 
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 02/2021


if ~exist('idx_without_NaN', 'var') || isempty(idx_without_NaN)
    idx_without_NaN = 2:length(X);
end

if ~exist('Y', 'var')
    Y = X;
end

% Initialize the variable count
lA = length(A);
ncontexts = length(contexts);
lw = cellfun(@(x) length(x), contexts);
height = max(lw);

T = {};     % past sequences
I = {};     % their positions in the sequence
Tidx = {};  % position in the list of contexts (0 if it is not a context)
tNwa = [];  % counts for the transitions to the next symbol

if ncontexts ~= 0
    for a = A
        [f, id, c, nwa] = find_context(a, lA, idx_without_NaN, X, contexts, ncontexts, height, Y);
        T = [T, f];
        I = [I, id];
        Tidx = [Tidx, c];
        tNwa = [tNwa; nwa];
    end
end

lc = length(Tidx);
count = cell(2, ncontexts + lc); %preallocate for speed
count(1,:) = {0};
Nwa = zeros(ncontexts + lc, lA);

% put the results in the structure of count
additional = 0;
for i = 1 : lc
    if Tidx{i} > 0
        count{1, Tidx{i}} = length(I{i}); 
        count{2, Tidx{i}} = I{i} + lw(Tidx{i}) - 1;
        try
        Nwa(Tidx{i},:) = tNwa(i,:);
        catch
            aa=1;
        end
    else
        % get the positions where the non-context appears
        positions = I{i} + length(T{i}) - 1;
        % add this non-context only if it appears at some position where is possible to find a context
        if any(positions >= height)    
            additional = additional + 1;
            count{1, ncontexts+additional} = length(I{i});
            count{2, ncontexts+additional} = positions;
            Nwa(ncontexts+additional,:) = tNwa(i,:); 
        end
    end
end
%shrink memory
count(:, ncontexts+additional+1:end) = [];
Nwa(ncontexts+additional+1:end,:) = [];

end