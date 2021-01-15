function [c, idx] = contextfunction(past, contexts)
%CONTEXTFUNCTION return the context associated to a given sequence (past sequence)
% 
% Inputs
%
% 	past      : sequence
% 	contexts  : set of contexts (context tree)
% 
% Outputs
%
% 	c         : context associated to the sequence [past]
% 	idx       : index (position) of the context [c] in [contexts]. If do
%               not exist, idx=-1

%
% Usage
%			ctxs = {0, [0 1], [1 1], 2};
%			past = [0 1 1 0 1 2 0 1];
%			[c, idx] = contextfunction(past, ctxs);
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 05/2020


n_contexts = length(contexts);
lpast = length(past);
idx = -1;
c = [];

% i = 1;
% li = length(contexts{i});
% while (i <= n_contexts)&&...
%         ( (li > lpast) || (~isequal(contexts{i},past(end-li+1 : end))) )
%     i = i + 1;
% end
% 
% if (i <= n_contexts)
%     c = contexts{i};
%     idx = i;
% end

found = false;
i = 0;
while (~found)&&(i < n_contexts)
    i = i + 1;
    li = length(contexts{i});
    found = (li <= lpast) && (isequal(contexts{i},past(end-li+1 : end)));
end

if found
    c = contexts{i};
    idx = i;
end

end