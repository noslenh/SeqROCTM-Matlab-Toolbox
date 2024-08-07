function [past, M, Mc, iT, past2ctx] = contextTree_to_FiniteMarkov(contexts, P, A)
%CONTEXTTREE_TO_FINITEMARKOV Compute the representation of a context tree model as a k-order Markov process.
%                            (k refers to the height of the tree)
% Inputs
%
%   contexts	: contexts of the context tree model
%   P       	: transition probabilities (each row is the distribution of a context)
%   A       	: alphabet
%
% Outputs
%
%   past    	: past sequences that the context tree model can generate
%   M       	: transition probabilities from one past to another 
%   Mc      	: transition probabilities of a symbol given a past
%   iT      	: matrix with the index of past that is formed given a past and a symbol
%   past2ctx	: array with the index of the context that corresponds to each past
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 08/2020

n_symbols = size(A, 2);
n_contexts = size(contexts, 2);

% get the contexts codified as indexes of symbols in the alphabet
for i = 1 : n_contexts
    [~, icontexts{i}] = ismember(contexts{i}, A);
end

% get the order of the Markov process (height of the tree)
L = max(cellfun(@(x) length(x), contexts, 'uniformoutput', 1));

% get all possible sequences of length L that can be formed with symbols in A
[all_past, I] = permn(A, L);     
n_past = size(I,1);

% check which pasts can be generated by the model (maybe there is a more
% efficient way of doing this)
existing_past = zeros(n_past, 1);
for i = 1 : n_past
    existing_past(i) = check_past_exist(I(i, :), icontexts, P);
end

M = zeros(n_past, n_past);       % transition prob. from past to past 
Mc = zeros(n_past, n_symbols);   % transition prob. from past to symbol
iT = zeros(n_past, n_symbols);   % index of past = past + symbol
past2ctx = zeros(n_past, 1);	 % index of the context associated to each past

% looks for the possible transitions for each past
for p = 1 : n_past 
    if existing_past(p)
        l_p = I(p,:);
        % context associated to l_p
        [~, iwp] = contextfunction(l_p, icontexts);
        [~, ipastT] = past_with_transitions(l_p, I, n_symbols);
        M(p, ipastT) = P(iwp, :);
        Mc(p, :) = P(iwp, :);
        % matrix with indexes of new_past = past + symbol
        iT(p, :) = ipastT;
		%
		past2ctx(p) = iwp;
    end
end

% delete from M the non_existing pasts
M(~existing_past,:) = [];
M(:,~existing_past) = [];
Mc(~existing_past,:) = [];         % trans. prob p(symbol|past)
past2ctx(~existing_past) = [];

past = all_past(logical(existing_past), :); % possible pasts 

cs = cumsum(~existing_past);
idp = (1:n_past)';

cs(~existing_past) = [];
idp(~existing_past) = [];

zz = zeros(n_past, 1);
zz(idp) = idp - cs; 

iT(~existing_past, :) = [];
iT = zz(iT);

end