function [HN, hn] = EntropyN(est_P, est_transP, iT, order)
% ENTROPYN computes the entropy of the distribution of subsequences of
%          length 1,...,N that are generated by a context tree model
%
% Inputs
% 
%  est_P         : cell array with estimates of probabilities of
%                  occurrence of all possible sequences of length 1,2,...,height. 
%  est_transP    : matrix with estimates of transition probabilities
%                  associated to sequences of length 'height'
%  iT            : cell array of matrices with the index of the past that is formed given a
%                  past and a symbol
%  order         : order of the Markov process (height of the context tree)
%
% Outputs
%
%  HN            : vector with the entropy of distribution of sequences of length 1,...,N
%  hn            : vector with block entropy
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 08/2020

N = 10;             % max length of sequence to be considered (the computational cost
                    % increase with this number) 

EN = zeros(N,1);

nsymbols = size(iT,2);

% compute the entropy of sequences of length 1,2,...,k using the
% probabilities
HN = zeros(N,1);
for k = 1 : order
    probs = est_P{k};
    probs(probs == 0) = [];
    HN(k) = -1 * sum(probs .* log2(probs));
end


k = order;
npastk = size(est_P{k}, 1);

for p = 1 : npastk
    for i = 1 : nsymbols
        EN = contribute(k + 1, est_P{k}(p), 0, p, i, est_transP, iT{k}, EN, N);
    end
end

HN(k+1: end) = HN(k) + EN(k+1:end);

% Assuming that H_0 = 0, we have
HN = [0; HN];
hn = diff(HN);

end

function EN = contribute(n, pn, slogn, lpast, s, transP, iT, EN, N)
%CONTRIBUTE Recursive function to assist the computation in EntropyN
%
% Inputs
%
%  n      : level of the entropy Hn this term is contributing to
%  pn     : p(s1,...,sn); s1,...,sn is the past from which this function was called
%  slogn  : sum_{i = k+1}^n log(p(s_i|s_{i-1}^{i-k})
%  lpast  : index of the past s1,...,sn
%  s      : index of the observed symbol
%  transP : matrix with the transition probabilities p(s|s_1,...,s_k)
%  iT     : matrix with the indexes of the past formed from a past and a symbol
%  EN     : 
%  N      : maximum length of subsequences  
%
% Outputs
%  EN     :   
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 08/2020

nsymbols = size(iT,2);
tp = transP(lpast,s); 

if (tp > 0)  % if this prob is zero, the result is zero
    % update entropy
    pn = pn * tp;
    slogn = slogn + log2(tp);

    EN(n) = EN(n) - pn*slogn;
    
    if (n < N) % if n = N, we do not need to call for the next level
        % call the following level
        for i = 1 : nsymbols
            EN = contribute(n + 1, pn, slogn, iT(lpast, s), i, transP, iT, EN, N);
        end
    end
end
end
