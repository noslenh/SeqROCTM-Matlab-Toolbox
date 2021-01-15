function [ctxs, index] = is_suffix_of_contexts(ss, contexts)
%IS_SUFFIX_OF_CONTEXTS Get all the contexts from which the sequence 'ss' is
%                      a suffix
% Inputs
%
% 	ss        : sequence
% 	contexts  : set of contexts
%
% Output
%
% 	ctxs      : contexts from which ss is a suffix
% 	index     : indexes of such contexts in the array contexts
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 05/2020

ncontexts = length(contexts);
ll = length(ss);

% allocate memory
index = -1*ones(ncontexts,1);
ctxs = cell(ncontexts,1);

nn = 0;
for i = 1 : ncontexts 
    if (length(contexts{i}) >= ll)&&(isequal(contexts{i}(end-ll+1 : end), ss)) %(sum(contexts{i}(end-ll+1 : end) == ss) == ll)
        nn = nn + 1;
        ctxs{nn} = contexts{i};
        index(nn) = i;
    end
end

% shrink the memory
index(nn+1:end) = [];
ctxs(nn+1:end) = [];