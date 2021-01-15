function [idx_opt_model, R] = modeltunning_risk(criteria, param_set, risk, B, n, tau0, P0, A)
%MODELTUNNING_RISK 
%
% Inputs
%
% criteria  : algorithm to estimate the context tree model
% param_set : set of values of the parameter to be tuned
% risk      : risk function 
% B         : number of bootstrap samples
% n         : length of the bootstrap samples 
% tau0      : context tree used to generate the bootstrap samples
% P0        : transition probabilities associated to the leaves of the tree
% A         : alphabet
%
% Outputs
%
% idx_opt_model : index of the optimal model
% R             : risk values corresponding to each parameter value 

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 04/2020


r = risk;

% height of the tree
max_height = max(cellfun(@(x) length(x), tau0, 'UniformOutput', 1));

% % initialize matrices
% SX = zeros(B,n);
% SY = zeros(B,n+1);

% Generate B bootstrap sample
%disp('Generating bootstrap samples...');
% SX = generatesampleCTM_fast(tau0, P0, A, n, B);
SY = generatesampleCTM_fast(tau0, P0, A, n+1, B);
SX = SY(:,1:n);

% store the index of the bootstrap sample with some incompatibility
bad_bootstrap_sample = [];

np = length(param_set);
L = zeros(B, np);
for k = 1 : np
%     disp(['Computing bootstrapped risk for parameter value ' num2str(param_set(k)) '...']);
    % estimate the model using param_set(k) for each bootstrap sample
    for b = 1 : B
        % estimate the model
        [ck, Pk] = estimate_contexttree(SX(b,:), A, max_height, criteria, param_set(k));
        % predict
        Yhat = predictor_delta_loss(SY(b,1:end-1), ck, Pk', A);
        
        if Yhat == -1
            % get the index of this bootstrap sample
            bad_bootstrap_sample = [bad_bootstrap_sample; b];
            disp(['Warning: The estimated model using Bootstrap sample ' int2str(b)...
                    ' was incompatible with the bootstrap sample generated for prediction.']);
        else
            % compute the loss 
            L(b,k) = delta_loss(Yhat, SY(b,end));
        end
    end
end
% bootstrap approximation of the Risk for each model
L(bad_bootstrap_sample,:) = []; %delete the information for the bootstrap samples with problem
R = sum(L)/(B-length(bad_bootstrap_sample));

% index of the best model
[~, tmp_idx] = min(R(end:-1:1));
idx_opt_model = np - tmp_idx + 1; 
end

function Xhat = predictor_delta_loss(X, contexts, P, A)
%PREDICTOR_DELTA_LOSS optimal predictor for the zero-one loss fuction
 
 if isempty(contexts)
     % i.i.d model, so return the most probable symbol
     [~, idxs] = max(P);
     Xhat = A(idxs);
 else
     % get the index of the context associated to the past
     [~, idxc] = contextfunction(X, contexts);
     if idxc == -1
        % the context was not found, so the model is not compatible
        % with the sample
        Xhat = -1;
     else
        % get the mode of the distribution associated to that context
        [~, idxs] = max(P(idxc,:));
        Xhat = A(idxs);
     end
 end
end

function L = delta_loss(X, Xhat)
    L = (X ~= Xhat);
end
