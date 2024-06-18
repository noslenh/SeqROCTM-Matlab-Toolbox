function logL = treeloglikelihood(X, tree, alphabet, missing)
%TREELOGLIKELIHOOD  Compute the likelihood of a context tree for the data X
% 
% Inputs
%
%   X           : sequence of inputs
%   tree        : context tree
%   alphabet    : alphabet 
%   missing     : 0 (false), 1 (true) or an array with the indexes without
%                   missing values              
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
    N = histc(X, alphabet);
    ind = N > 0;
    logL = sum(N(ind) .* (log(N(ind)) - log(nsample)));
else
    % Before counting, analyze if there are missing values
    if ~exist('missing', 'var') || isequal(missing, 0)
        % do the counting using all the sequence (as usual)
        [Count, N] = countctx(tree, X, alphabet);
    else
        if missing == 1
            % calculate the positions without NaN (excluding the first position)
            idx_without_NaN = find(~isnan(X(2:end)));
            % sum 1 to update the indexes
            idx_without_NaN = idx_without_NaN + 1;
        else
            % get the positions without NaN
            idx_without_NaN = missing;
        end
        [Count, N] = countctx(tree, X, alphabet, idx_without_NaN);
    end

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