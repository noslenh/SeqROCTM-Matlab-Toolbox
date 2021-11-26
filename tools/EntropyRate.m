function [H, Hmu, mu] = EntropyRate(S)
%ENTROPYRATE Compute the entropy rate of a finite order Markov process
%            with transition probabilities between the states specified by
%            the matrix S
%
% Inputs
%
% S    : Matrix specifying the transition probabilities between the
%        states of the Markovian system
%
% Outputs
%
% H    : entropy rate of the process
% Hmu  : entropy of the stationary distribution
% mu   : probability of occurrence of each state (stationary distribution)


%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 08/2020
    
no_states = size(S,1);

% compute the stationary distribution
[V, D] = eig(S');
index = find((abs(diag(D) - 1) < 10^-10));
mu = V(:,index)/sum(V(:,index));


H = 0;
Hmu = 0;
for i = 1 : no_states
    % entropy rate
    T = 0;
    for j = 1 : no_states
        if S(i,j) ~= 0, T = T + S(i,j)*log2(S(i,j)); end
    end
    H = H + mu(i)*(-T);
    
    % entropy of the stationary distribution
    if mu(i) ~= 0, Hmu = Hmu - mu(i)*log2(mu(i)); end
end
    
end