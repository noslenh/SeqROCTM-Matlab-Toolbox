function [idx_opt_model, R] = tuning_risk2(param_set, Xbootsamples, Ybootsamples, A, options)
%MODELTUNNING_RISK2 Model selection using a risk function for a SeqROCTM
%
% Inputs
%
%   param_set      : set of values of the parameter to be tuned
%   Xbootsamples   : bootstrap samples for the sequence X
%   Ybootsamples   : bootstrap samples for the sequence Y
%   A              : alphabet
%   options        : structure with the values to be passed to the
%                       estimation functions
%
% Outputs
%
% idx_opt_model : index of the optimal parameter value
% R             : risk values corresponding to the parameter values 

%   References:
%      [1] N. Hernández et al., arXiv:2009.06371, (2021).  
%      [2] P. Buhlmann et al., Ann. Inst. Statist. Math, 52, 1, 287-315 (2000)
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 01/2021

lA = length(A);
B = size(Xbootsamples,1);

% index of the bootstrap sample with some incompatibility
bad_bootstrap_sample = [];

np = length(param_set);
L = zeros(B, np);

if any(strcmpi(options.EstimationMethod, {'context_empD', 'context_cL'}))
    %%% Context Algorithm
    for b = 1 : B
        
        Xb = Xbootsamples(b,1:end-1);
        Yb = Ybootsamples(b,1:end-1);
        
        % some computations to speed-up
        % compute the complete tree and the TEST structure only once (for speed-up)
        [T, I] = completetree(Xb, options.MaxTreeHeight, A);
        TEST = getTESTstructure(T, I, lA, Yb);
        ct_inf = {T, I};
        %BIC-info
        precomputed_stats = [];
        
        % for each value param_set(k) of the parameter and each bootstrap
        % sample
        for k = 1 : np
            %estimate the context tree model
            [ck, Pk] = estimate_discreteSeqROCTM(Xb, Yb, A, ...
                'MaxTreeHeight', options.MaxTreeHeight,             ...
                'EstimationMethod', options.EstimationMethod,       ...
                'ParameterValue', param_set(k),                     ...
                'CtxCompleteTree', ct_inf,                          ...
                'CtxTestStructure', TEST,                           ...
                'BicDegreeOfFreedom', options.BicDegreeOfFreedom,   ...
                'BicPrecomputedStats', precomputed_stats);
            
            % predict
            Yhat = predictor_delta_loss(Xb, ck, Pk, A);
            
            % compute the loss
            if Yhat == -1
                % get the index of this bootstrap sample (this can happens if
                % the length of the sequences is small)
                bad_bootstrap_sample = [bad_bootstrap_sample; b];
                disp(['Warning: The estimated model using Bootstrap sample ' int2str(b)...
                    ' was incompatible with the bootstrap sample generated for prediction.']);
            else
                % loss
                L(b,k) = delta_loss(Yhat, Xbootsamples(b,end));
            end
        end
    end
elseif strcmpi(options.EstimationMethod, 'bic')
    %%% BIC Algorithm
    for b = 1 : B
        
        Xb = Xbootsamples(b,1:end-1);
        Yb = Ybootsamples(b,1:end-1);
        
        % some computations to speed-up
        % compute the statistics Nw and Nwa used in BIC only once (for speed-up)
        df1 = ~strcmpi(options.BicDegreeOfFreedom,'fix');
        [~, ~, ~, outps] = bic_WCT(Xb, A, options.MaxTreeHeight, param_set(1), df1, 0, Yb);
        precomputed_stats{1} = outps.stats(:,5:6+lA-1);
        precomputed_stats{2} = outps.nonExistingNodes;
        precomputed_stats{3} = outps.XlengthWithoutNaN;
        %ctx-info
        TEST = -1;
        ct_inf = -1;
        
        % for each value param_set(k) of the parameter and each bootstrap
        % sample
        for k = 1 : np
            %estimate the context tree model
            [ck, Pk] = estimate_discreteSeqROCTM(Xb, Yb, A, ...
                'MaxTreeHeight', options.MaxTreeHeight,             ...
                'EstimationMethod', options.EstimationMethod,       ...
                'ParameterValue', param_set(k),                     ...
                'CtxCompleteTree', ct_inf,                          ...
                'CtxTestStructure', TEST,                           ...
                'BicDegreeOfFreedom', options.BicDegreeOfFreedom,   ...
                'BicPrecomputedStats', precomputed_stats);
            
            % predict
            Yhat = predictor_delta_loss(Xb, ck, Pk, A);
            
            % compute the loss
            if Yhat == -1
                % get the index of this bootstrap sample (this can happens if
                % the length of the sequences is small)
                bad_bootstrap_sample = [bad_bootstrap_sample; b];
                disp(['Warning: The estimated model using Bootstrap sample ' int2str(b)...
                    ' was incompatible with the bootstrap sample generated for prediction.']);
            else
                % loss
                L(b,k) = delta_loss(Yhat, Xbootsamples(b,end));
            end
        end
    end
else
    error('The estimation method %s is not implemented to be tunned', options.EstimationMethod);
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
