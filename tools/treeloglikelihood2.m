function logL = treeloglikelihood2(X, Y, tree, alphabet)
%TREELOGLIKELIHOOD  Compute the likelihood of a context tree for the SeqROCTM (X,Y)
% Inputs
%
%   X           : sequence of inputs
%   Y           : sequence of responses (optional)
%   tree        : context tree
%   alphabet    : alphabet 
%
% Outputs
%
%   logL        : log likelihood
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 01/2021

ncontexts = length(tree);
nsymbols = length(alphabet);
nsample = length(X);

if isempty(tree)  
    % return the entropy since the sequence is iid
    N = histc(Y, alphabet);
    ind = N > 0;
    logL = sum(N(ind) .* (log(N(ind)) - log(nsample)));
else
    [Count, N] = countctx(tree, X, alphabet, [], Y);
    if size(Count,2) > ncontexts %there are pasts in the sequence that have no context associated,
        logL = -inf;             %so the likelihood of the model is zero   
    else
        % compute the log-likelihood
        ss = cell2mat(Count(1,:));
                
        % elements different from zero
        ind = N > 0;
        B = repmat(ss', 1, nsymbols);

        % log-likelihood
        logL = sum( N(ind) .* ( log(N(ind)) - log(B(ind)) ) );
    end
end