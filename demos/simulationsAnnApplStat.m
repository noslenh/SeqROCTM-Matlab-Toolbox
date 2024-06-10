% This script reproduces the simulations done in the article [A. Galves
% et. al., Ann. Appl. Stat., Volume 6, Number 1 (2012), 186-209] using the
% Matlab SeqROCTM package.

% In the simulations, the symbol 1 was used as a renewal point to generate
% the bootstrap samples. 

% Here we include the routines to do also the simulations with the Context
% Algorithm.

%%% Model specification %%%

% alphabet
A = [0,1];

% contexts
contexts  = {1,  [1 0], [1 0 0], [0 0 0]};

% % family of distributions - Model 1
% P = [1, 0; 0.3, 0.7; 0.2, 0.8; 0.25, 0.75];  

% family of distributions - Model 2
P = [1, 0; 0.2, 0.8; 0.3, 0.7; 0.4, 0.6]; 

% renewal point specified by the user
th_renwpoint = 1;


%%% Parameters value %%%%

n 			= 10000;			% length of the stochastic sequence
Repetitions = 100;				% number of times the procedure is repeated
B 			= 200;				% number of bootstrap samples
n1 			= floor(0.3*n); 	% proportion of the size of the sample corresponding to the size of the smaller re-sample.
n2 			= floor(0.9*n);		% proportion of the size of the sample corresponding to the size of the larger re-sample.
alpha 		= 0.01;            	% alpha level to use on the t-test
max_height 	= 6;				% height of the complete tree
c_min 		= 0;				% minimum value of the BIC constant
c_max 		= 1000;				% maximum value of the BIC constant
c_max_ctx   = c_max*log(n);

% % fix the seed if you want to control random generations
% rng(200);

%%% Simulations %%%%

% number of contexts
ncontexts = numel(contexts);

% initialization 
inside_champions_bic = 0;
inside_champions_ctx = 0;
true_model_bic = 0;
true_model_ctx = 0;

% for each repetition
for r = 1 : Repetitions
    
    % generate a sequence of length n
%     X = generatesampleCTM(contexts, P, A, n);
    X = model2_10000(r,:);
    
    disp(['Processing sample ' num2str(r) ' ...']);
  
    % tune the context tree model using bic
    [optmodel, ~, resultsb] = tune_contextTreeModel(X, A, 'MaxTreeHeight', max_height,      ...
                                                          'ParameterLowerBound', c_min,     ...
                                                          'ParameterUpperBound', c_max,     ...
                                                          'BootRenewalPoint', th_renwpoint, ...
                                                          'BootNsamples', B,                ...
                                                          'n1', n1,                         ...
                                                          'n2', n2                          ...
                                                    );
    

    % check if the true model is within the Champion Trees and if it was
    % chosen as optimal
    if isequalCT(contexts, optmodel)
        true_model_bic = true_model_bic + 1;
        inside_champions_bic = inside_champions_bic + 1;
    else
        nl = cellfun(@(x) length(x), resultsb.champions);
        idx = find(nl == ncontexts);
        if (~isempty(idx))&&(isequalCT(contexts, resultsb.champions{idx}))
            inside_champions_bic = inside_champions_bic +1;
        end
    end
    
    % tune the context tree model using Context Algorithm
    [optmodel, ~, resultsc] = tune_contextTreeModel(X, A, 'MaxTreeHeight', max_height,      ...
                                                          'ParameterLowerBound', c_min,     ...
                                                          'ParameterUpperBound', c_max_ctx, ...
                                                          'BootRenewalPoint', th_renwpoint, ...
                                                          'BootNsamples', B,                ...
                                                          'n1', n1,                         ...
                                                          'n2', n2,                         ...
                                                          'EstimationMethod', 'context'     ...
                                                    );
                                                
    % check if the true model is within the Champion Trees and if it was
    % chosen as optimal
    if isequalCT(contexts, optmodel)
        true_model_ctx = true_model_ctx + 1;
        inside_champions_ctx = inside_champions_ctx + 1;
    else
        nl = cellfun(@(x) length(x), resultsc.champions);
        idx = find(nl == ncontexts);
        if (~isempty(idx))&&(isequalCT(contexts, resultsc.champions{idx}))
            inside_champions_ctx = inside_champions_ctx +1;
        end
    end
end

disp(['BIC: True model inside the Champion Trees: ' num2str(inside_champions_bic)]);
disp(['BIC: True model chosen: ' num2str(true_model_bic)]);
disp(['Context Algorithm: True model inside the Champion Trees: ' num2str(inside_champions_ctx)]);
disp(['Context Algorithm: True model chosen: ' num2str(true_model_ctx)]);
