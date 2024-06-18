function [pasts, indexes] = past_with_transitions(past, all_past, N)
%PAST_WITH_TRANSITIONS Return from a list of pasts the ones that can be generated from a given past.
%                      This function gives all the pasts (and its indexes)
%                      from all_past that can be generated from past.
% Inputs
%  
%  past 	: past sequence	
%  all_past : all past sequences
%  N		:
%
% Outputs
%
%  pasts 	: all pasts with transition from past
%  indexes 	: indexes of pasts in all_past
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 08/2020

L = length(past);
sub_past = past(2:end);

% looks for all the past that begin with sub_past (for this we made a computation
% that gets directly the indexes of those pasts with no need to iterate all_past)
index = 0;
for i = 1 : L - 1
    index = index + (sub_past(i) - 1) * perm_with_repetition_ordermatter(N, L-i);
end
indexes = index + (1: N);
pasts = all_past(indexes,:);
end

function n = perm_with_repetition_ordermatter(N, k)
    n = N^k;
end

