% Using the scenarios simulated in A. Galves et. al., Ann. Appl. Stat., 6,
% 1, 186-209 (2012) this script execute the different tuning procedures
% implemented in the SeqROCTM toolbox.

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

n 			= 10000; 			% length of the stochastic sequence
Repetitions = 100;				% number of times the procedure is repeated
B 			= 200;				% number of bootstrap samples
n1 			= floor(0.3*n); 	% proportion of the size of the sample corresponding to the size of the smaller resample.
n2 			= floor(0.9*n);		% proportion of the size of the sample corresponding to the size of the larger resample.
alpha 		= 0.01;            	% alpha level to use on the t-test
max_height 	= 6;				% height of the complete tree
c_min 		= 0;				% minimum value of the BIC constant
c_max 		= 1000;				% maximum value of the BIC constant
c_max_ctx   = c_max*log(n);

options.EstimationMethod = 'bic';
options.MaxTreeHeight = max_height;
options.BicDegreeOfFreedom = 'fix';
options.BicMissing = 0;

% % fix the seed if you want to control random generations
% rng(200);

%%% Simulations %%%%

% number of contexts
ncontexts = numel(contexts);

% indexes of the model chosen (for each sample and each of the 18 tuning
% configurations)
idx_optimalModel = zeros(100, 18);
champions_bic_fix = cell(100,1);
champions_bic_var = cell(100,1);
champions_ctx = cell(100,1);

% for each repetition
for r = 1 : Repetitions
    
    % generate a sequence of length n
    X = generatesampleCTM(contexts, P, A, n);
    
    disp(['Processing sample ' num2str(r) ' ...']);
  
    %%% CHAMPION TREES -> BIC - df=fix
    [champions, Ps, ~, prmvalues] = estimate_championTrees(X, A, 'MaxTreeHeight', max_height, ...
                                                        'ParameterLowerBound', 0, 'ParameterUpperBound', c_max);
    
    options.EstimationMethod = 'bic';
    options.BicDegreeOfFreedom = 'fix';
    
    % Bootstrap samples
    % the bootstrap samples are generate with length n+1, and shorten
    % accordingly
    
    % block - given renewal point
    b_blocks_g = bootstrap_blocks(X, th_renwpoint, n+1, B);
    
    % blocks - computed renewal point
    renewal_point = tree_renewalpoint(champions{1}, Ps{1}, A, X);
    b_blocks_c = bootstrap_blocks(X, renewal_point, n+1, B);
    
    % parametric
    b_param = generatesampleCTM_fast(champions{1}, Ps{1}, A, n+1, B);
    
    %BIC - df=fix - smc - block (given renw. point)
    [~, opt] = tuning_SMC(champions, A, n1, n2, alpha, b_blocks_g(:,1:n2), 0);
    idx_optimalModel(r,1) = opt;
    
    %BIC - df=fix - smc - block (computed renw. point)
    [~, opt] = tuning_SMC(champions, A, n1, n2, alpha, b_blocks_c(:,1:n2), 0);
    idx_optimalModel(r,2) = opt;
    
    %BIC - df=fix - smc - parametric
    [~, opt] = tuning_SMC(champions, A, n1, n2, alpha, b_param(:,1:n2), 0);
    idx_optimalModel(r,3) = opt;
    
    %BIC - df=fix - risk - block (given renw. point)
    opt = tuning_risk(prmvalues, b_blocks_g, A, options);
    idx_optimalModel(r,4) = opt;
    
    %BIC - df:fix - risk - block (computed renw. point)
    opt = tuning_risk(prmvalues, b_blocks_c, A, options);
    idx_optimalModel(r,5) = opt;
    
    %BIC - df:fix - risk - parametric
    opt = tuning_risk(prmvalues, b_param, A, options);
    idx_optimalModel(r,6) = opt;
    
    champions_bic_fix{r} = champions;    
    
    %%% CHAMPION TREES -> BIC - df=variable
    [champions, Ps, ~, prmvalues] = estimate_championTrees(X, A, 'MaxTreeHeight', max_height, ...
                            'ParameterLowerBound', 0, 'ParameterUpperBound', c_max, ...
                            'BicDegreeOfFreedom', 'variable');
    
    options.EstimationMethod = 'bic';
    options.BicDegreeOfFreedom = 'variable';
                        
    % Bootstrap samples
    % blocks - computed renewal point
    renewal_point = tree_renewalpoint(champions{1}, Ps{1}, A, X);
    b_blocks_c = bootstrap_blocks(X, renewal_point, n+1, B);
    
    % parametric
    b_param = generatesampleCTM_fast(champions{1}, Ps{1}, A, n+1, B);
    
    %BIC - df:variable - smc - block (given renw. point)
    [~, opt] = tuning_SMC(champions, A, n1, n2, alpha, b_blocks_g(:,1:n2), 0);
    idx_optimalModel(r,7) = opt;
    
    %BIC - df:variable - smc - block (computed renw. point)
    [~, opt] = tuning_SMC(champions, A, n1, n2, alpha, b_blocks_c(:,1:n2), 0);
    idx_optimalModel(r,8) = opt;
    
    %BIC - df:variable - smc - parametric
    [~, opt] = tuning_SMC(champions, A, n1, n2, alpha, b_param(:,1:n2), 0);
    idx_optimalModel(r,9) = opt;
    
    %BIC - df:variable - risk - block (given renw. point)
    opt = tuning_risk(prmvalues, b_blocks_g, A, options);
    idx_optimalModel(r,10) = opt;
    
    %BIC - df:variable - risk - block (computed renw. point)
    opt = tuning_risk(prmvalues, b_blocks_c, A, options);
    idx_optimalModel(r,11) = opt;
    
    %BIC - df:variable - risk - parametric
    opt = tuning_risk(prmvalues, b_param, A, options);
    idx_optimalModel(r,12) = opt;
    
    champions_bic_var{r} = champions; 
    
    
    % CHAMPION TREES -> Context Algorithm
    [champions, Ps, ~, prmvalues] = estimate_championTrees(X, A, 'MaxTreeHeight', max_height, ...
                                            'ParameterLowerBound', 0, 'ParameterUpperBound', c_max_ctx, ...
                                            'BicDegreeOfFreedom', 'variable', ...
                                            'EstimationMethod', 'context');
                                        
    options.EstimationMethod = 'context';
                    
    % blocks - computed renewal point
    renewal_point = tree_renewalpoint(champions{1}, Ps{1}, A, X);
    b_blocks_c = bootstrap_blocks(X, renewal_point, n+1, B);
    
    % parametric
    b_param = generatesampleCTM_fast(champions{1}, Ps{1}, A, n+1, B);
    
    %CTX - smc - block (given renw. point)
    [~, opt] = tuning_SMC(champions, A, n1, n2, alpha, b_blocks_g(:,1:n2), 0);
    idx_optimalModel(r,13) = opt;
    
    %CTX - smc - block (computed renw. point)
    [~, opt] = tuning_SMC(champions, A, n1, n2, alpha, b_blocks_c(:,1:n2), 0);
    idx_optimalModel(r,14) = opt;
    
    %CTX - smc - parametric
    [~, opt] = tuning_SMC(champions, A, n1, n2, alpha, b_param(:,1:n2), 0);
    idx_optimalModel(r,15) = opt;
    
    %CTX - risk - block (given renw. point)
    opt = tuning_risk(prmvalues, b_blocks_g, A, options);
    idx_optimalModel(r,16) = opt;
    
    %CTX - risk - block (computed renw. point)
    opt = tuning_risk(prmvalues, b_blocks_c, A, options);
    idx_optimalModel(r,17) = opt;
    
    %CTX - risk - parametric
    opt = tuning_risk(prmvalues, b_param, A, options);
    idx_optimalModel(r,18) = opt;
    
    champions_ctx{r} = champions;
 
end

%%%% Summarization of the results
TM_inside_champions = zeros(3,1); % number of times the true model was inside the Champion Trees
TM_choosen = zeros(100,18);

%BIC-fix
for r = 1 : Repetitions
    nl = cellfun(@(x) length(x), champions_bic_fix{r});
    idx = find(nl == 4);
    if (~isempty(idx))&&(isequalCT(contexts, champions_bic_fix{r}{idx}))
        TM_inside_champions(1) = TM_inside_champions(1) + 1;
        TM_choosen(r,1:6) = (idx_optimalModel(r,1:6) == idx);
    end     
end

%BIC-variable
for r = 1 : Repetitions
    nl = cellfun(@(x) length(x), champions_bic_var{r});
    idx = find(nl == 4);
    if (~isempty(idx))&&(isequalCT(contexts, champions_bic_var{r}{idx}))
        TM_inside_champions(2) = TM_inside_champions(2) + 1;
        TM_choosen(r,7:12) = (idx_optimalModel(r,7:12) == idx);
    end     
end

%Context
for r = 1 : Repetitions
    nl = cellfun(@(x) length(x), champions_ctx{r});
    idx = find(nl == 4);
    if (~isempty(idx))&&(isequalCT(contexts, champions_ctx{r}{idx}))
        TM_inside_champions(3) = TM_inside_champions(3) + 1;
        TM_choosen(r,13:18) = (idx_optimalModel(r,13:18) == idx);
    end     
end

% Number of time each method choose the true model
nselected = sum(TM_choosen);
