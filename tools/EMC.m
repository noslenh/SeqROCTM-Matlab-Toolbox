function compl = EMC(contexts, P, A)
%EMC Compute the Effective Measure Complexity of a context tree model (see
%    Grassberger, International Journal of Theoretical Physics 25, 1986) 
%
% Inputs
%
%  contexts : set of contexts
%  P        : transition probabilities associated to the contexts
%  A        : alphabet
%
% Output
%
%  compl    : Effective measure complexity
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 08/2020

order = max(cellfun(@(x) length(x), contexts, 'uniformoutput', 1));

% estimate some probabilities and transition probabilities
[est_P, est_transP, iT] = empprobsubsequences(contexts, P, A);

% get the block entropy
[~, hn] = EntropyN(est_P, est_transP, iT, order);

% get the entropy rate
H = EntropyRateCT(contexts, P, A);

% compute the complexity EMC
compl = 0;
for i = 1 : order
    compl = compl + hn(i) - H;
end

end