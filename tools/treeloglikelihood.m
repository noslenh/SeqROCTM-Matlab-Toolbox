function logL = treeloglikelihood(tree, X, alphabet)
%TREELOGLIKELIHOOD  Compute the likelihood of a context tree for the data X
% Inputs
%
%   tree        : context tree
%   X           : sequence of symbols taking values in the alphabet
%   alphabet    : alphabet 
%
% Outputs
%
%   logL        : log likelihood
%
% Usage
%			A = [0,1,2];
%
%			ctxs = {0, [0 1], [1 1], 2}
%			P = [0,   1,   0; ...                  
%     			 0,   0.25,   0.75;
%     			 1,   0,   0;
%     			 1,   0,   0 ];
%
%			X = generatesampleCTM(ctxs, P, A, 100);  
%			logL = completetree(ctxs, X, A);
% 

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 05/2019

ncontexts = length(tree);
nsymbols = length(alphabet);
nsample = length(X);

if isempty(tree)  
    % return the entropy since the sequence is iid
    N = histc(X, alphabet);
    ind = N > 0;
    logL = sum(N(ind) .* (log(N(ind)) - log(nsample)));
else
    Count = countctx(tree, X, alphabet);
    if size(Count,2) > ncontexts %there are pasts in the sequence that have no context associated,
        logL = -inf;             %so the likelihood of the model is zero   
    else
        % log-likelihood
        N = zeros(ncontexts, nsymbols);
        ss = zeros(ncontexts, 1);

        for i = 1 : ncontexts
            idx = Count{2,i};
            ss(i) = Count{1,i};
            for id = idx
                column = X(id + 1) + 1; %trick: use the fact that the alphabet is always [0,...,|A|-1]
                N(i,column) = N(i,column) + 1; 

            end 
        end
        % elements different from zero
        ind = N > 0;
        B = repmat(ss, 1, nsymbols);

        % log-likelihood
        logL = sum( N(ind) .* ( log(N(ind)) - log(B(ind)) ) );
    end
    
end