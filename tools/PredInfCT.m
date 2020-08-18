function [predinf, H] = PredInfCT(contexts, P, A)
%IPRED Compute Predictive Information of a context tree model (see Bialek
%      et.al, Neural Computation 13, 2409-2463, 2001)
%
% Inputs
%
%  contexts : set of contexts
%  P        : transition probabilities associated to the contexts
%  A        : alphabet
%
% Output
%
%  predinf  : Predictive information
%  H        : Entropy rate 
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 08/2020

% height of the context tree
height = max(cellfun(@(x) length(x), contexts, 'uniformoutput', 1));

% compute the entropy rate and the stationary distribution of the finite
% order Markov process corresponding to the context tree model
[H, ~, ~, Hmu_markovian] = EntropyRateCT(contexts, P, A);

% predictive information (this is an specific derivation for Markov process
% of finite order)
predinf = Hmu_markovian - height*H;

end
