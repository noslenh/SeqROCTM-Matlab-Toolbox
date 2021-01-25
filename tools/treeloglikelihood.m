function logL = treeloglikelihood(X, tree, alphabet, missing)
%TREELOGLIKELIHOOD  Compute the likelihood of a context tree for the data X
%                   or for the SeqROCTM (X,Y)
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
%			logL = completetree(ctxs, A, X);
% 

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
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
    % Before counting, analize if there are missing values
    if ~exist('missing', 'var') || isequal(missing, 0)
        % do the counting using all the sequence (as ussual)
        Count = countctx(tree, X, alphabet);
    else
        if missing == 1
            % calculate the positions without NaN (exluding the first position)
            idx_without_NaN = find(~isnan(X(2:end)));
            % sum 1 to update the indexes
            idx_without_NaN = idx_without_NaN + 1;
        else
            % get the positions without NaN
            idx_without_NaN = missing;
        end
        Count = countctx(tree, X, alphabet, idx_without_NaN);
    end

%     Count = countctx(tree, X, alphabet);
    if size(Count,2) > ncontexts %there are pasts in the sequence that have no context associated,
        logL = -inf;             %so the likelihood of the model is zero   
    else
        % compute the log-likelihood
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