function [past, indexes] = past_with_transitions(past, all_past, N)

L = length(past);
sub_past = past(2:end);

% looks for all the past that begin with sub_past (esto saca una cuenta que
% te ubica directamente en el íncide de esos pasados, sin necesidad de
% recorrer todos los pasados)
index = 0;
for i = 1 : L - 1
    index = index + (sub_past(i) - 1) * perm_with_repetition_ordermatter(N, L-i);
end
indexes = index + (1: N);
past = all_past(indexes,:);
end

function n = perm_with_repetition_ordermatter(N, k)
    n = N^k;
end

