function [H, Hmu, mu_ctx, Hmu_markovian] = EntropyRateCT(contexts, P, A)
%ENTROPYRATECT Compute the entropy rate of the stochastic process
%              represented by the context tree model. This is done by
%              representing the context tree model as a finite order Markov
%              chain and computing the entropy rate of such a Markovian
%              system
%
% Input
%
% contexts      : set of contexts
% P             : transition probabilities associated to the contexts
% A             : alphabet
%
% Output
%
% H             : entropy rate of the process
% Hmu           : entropy of the stationary distribution
% mu_ctx        : probability of occurrence of contexts (stationary distribution)
% Hmu_markovian : entropy of the stationary distribution of the Markovian system

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 08/2020

% context tree model to finite order Markov chain
[~, M, ~, ~, past2ctx] = contextTree_to_FiniteMarkov(contexts, P, A);

% compute the entropy rate 
[H, Hmu_markovian, mu] = EntropyRate(M);

% compute the stationary distribution of contexts
ns = size(M,1);
mu_ctx = zeros(length(contexts), 1);
for s = 1 : ns
    mu_ctx(past2ctx(s)) = mu_ctx(past2ctx(s)) + mu(s);
end

% entropy of the distribution mu_ctx
Hmu = sum(-mu_ctx.*log2(mu_ctx));

end
