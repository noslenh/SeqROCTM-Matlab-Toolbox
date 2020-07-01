% This script reproduces the simulations done in the article [A. Galves
% et. al., Ann. Appl. Stat., Volume 6, Number 1 (2012), 186-209] using the
% Matlab SeqROCTM package.

% In the simulations done in the article the symbol 1 was used as a renewal
% point to generate the bootstrap samples. Here we included in the
% simulations two other bootstrap strategies implemented in the package
% than can be useful in realistic situations (true model and renewal points
% unknown)

%%% Model specification %%%

% alphabet
A = [0,1];

% contexts
contexts  = {1,  [1 0], [1 0 0], [0 0 0]};

% family of distributions - Model 1
P = [1, 0; 0.3, 0.7; 0.2, 0.8; 0.25, 0.75];  

% % family of distributions - Model 2
% P = [1, 0; 0.2, 0.8; 0.3, 0.7; 0.4, 0.6]; 

% renewal point specified by the user
th_renwpoint = 1;


%%% Parameters value %%%%

Repetitions = 100;				% number of times the procedure is repeated
B 			= 200;				% number of bootstrap samples
n 			= 5000;				% length of the stochastic sequence
n1 			= floor(0.3*n); 	% proportion of the size of the sample corresponding to the size of the smaller resample.
n2 			= floor(0.9*n);		% proportion of the size of the sample corresponding to the size of the larger resample.
alpha 		= 0.01;            	% alpha level to use on the t-test
max_height 	= 6;				% height of the complete tree
c_min 		= 0;				% minimum value of the BIC constant
c_max 		= 500;				% maximum value of the BIC constant

% fix the seed if you want to control random generations
rng(200);

%%% Simulations %%%%

% number of contexts
ncontexts = numel(contexts);

% initialization 
mG1 = []; mG2 = []; mG3 = [];

% for each repetition
for r = 1 : Repetitions
    
    %generate a sequence of length n
    X = generatesampleCTM(contexts, P, A, n);
  
    %estimate the champion trees
    [trees, Ps, ML, cutoff] = estimate_championTrees(X, max_height, A, c_min, c_max);
    
% 	%%% if you want to plot the curve, uncomment this %%%
%     % plot the curve models vs. Likelihood
%      figure
%      nleaves = cellfun(@(x) size(x,2), trees);
%      plot(nleaves, ML, '*--b')
%      ylabel('log-likelihood');
%      xlabel('no. of contexts');

    % check if the true model is within the Champion Trees (this is only to speed-up the simulations:
	% if the Champion Trees does not contain the true model it is not necessary to carry out the procedure)
    nl = cellfun(@(x) length(x), trees);
    idx = find(nl == ncontexts);
    if (~isempty(idx))&&(isequalCT(contexts, trees{idx}))
        
        % we call the same procedure using different bootstrap strategies
        
		% bootstrap using an a priori known renewal point (this is the strategy used in the article)
		[~, id_G1] = modeltunning_championTrees(trees, A, n1, n2, alpha, B, 'blocks', X, th_renwpoint);
        
		% bootstrap finding a renewal point of the largest model in the Champion Trees 
        renewalpoint = tree_renewalpoint(trees{1}, Ps{1}', A, X);
        [~, id_G2] = modeltunning_championTrees(trees, A, n1, n2, alpha, B, 'blocks', X, renewalpoint);

		% parametric bootstrap using the largest model in the Champion Trees
        [~, id_G3] = modeltunning_championTrees(trees, A, n1, n2, alpha, B, 'parametric_ctm', trees{1}, Ps{1}');
		
        % check if the resulting context tree matches the true context tree
		mG1 = [mG1; isequalCT(trees{id_G1}, contexts)];
        mG2 = [mG2; isequalCT(trees{id_G2}, contexts)];
        mG3 = [mG3; isequalCT(trees{id_G3}, contexts)];
    end
end

disp(['The Champion Trees contained the true model: ' num2str(length(mG1)*100/Repetitions) ' %']);
disp(['Strategy 1: The true model was recovered ' num2str(sum(mG1)*100/Repetitions) ' %' ]);
disp(['Strategy 2: The true model was recovered ' num2str(sum(mG2)*100/Repetitions) ' %' ]);
disp(['Strategy 3: The true model was recovered ' num2str(sum(mG3)*100/Repetitions) ' %' ]);
