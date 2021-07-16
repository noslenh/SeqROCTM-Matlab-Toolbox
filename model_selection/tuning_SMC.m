function [opt_tree, idx] = tuning_SMC(championTrees, A, n1, n2, alpha, bootsamples_n2, missing)
%MODELTUNNING_SMC Context tree selection using the Smallest Maximizer Criterion 
%
% Inputs
%
%   championTrees   : set of Champion Trees (trees obtained for different
%                       values of the penalization constant in the BIC
%                       criteria)
%   A               : alphabet
%   n1              : proportion of the size of the sample corresponding to
%                       the size of the smaller resample.
%   n2              : proportion of the size of the sample corresponding to
%                       the size of the larger resample.
%   alpha           : significance level used in the t-test
%   bootsamples_n2  : bootstrap samples
%   missing         : 1 if treatment of Nan values is needed, 0 otherwise.
%
% Outputs
%
%  opt_tree         : optimal context tree
%  idx              : index of the optimal context tree in the set of
%                       Champion Trees
%

%   References:
%      [1] A. Galves et al., Ann. Appl. Stat., 6, 1, 186-209 (2012)


%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 01/2021

nTrees = length(championTrees);

B = size(bootsamples_n2,1);

if missing  %take into account that the data has missing values
    param_likhd_n1 = cell(B,1);
    param_likhd_n2 = cell(B,1);
    % compute the non_nan_indexes for the Bootstrap samples
    for b = 1 : B
        param_likhd_n1{b} = find(~isnan(bootsamples_n2(b,1:n1)));
        param_likhd_n2{b} = [param_likhd_n1{b}, n1+find(~isnan(bootsamples_n2(b,n1+1:end)))];
    end
else
    %initialize in such a way that always call the likelihood function
    %with missing in false
    param_likhd_n1(1:B,1) = {0};
    param_likhd_n2(1:B,1) = {0};
end

% compute the differences in likelihood for each pair of consecutive trees
% and all the bootstrap samples
diff_n1 = zeros(nTrees-1, B);
diff_n2 = zeros(nTrees-1, B);

% initialize L_current
L_current = zeros(B,2);
for b = 1 : B
    L_current(b,1) = treeloglikelihood(bootsamples_n2(b, 1:n1), championTrees{1}, A, param_likhd_n1{b});
    L_current(b,2) = treeloglikelihood(bootsamples_n2(b,:), championTrees{1}, A, param_likhd_n2{b});
end

for t = 1 : nTrees-1
    L_next = zeros(B,2); % store the log-likelihood of tree t+1 to speed-up
    for b = 1 : B
        %
        L_next(b,1) = treeloglikelihood(bootsamples_n2(b, 1:n1), championTrees{t+1}, A, param_likhd_n1{b});
        L_next(b,2) = treeloglikelihood(bootsamples_n2(b,:), championTrees{t+1}, A, param_likhd_n2{b});
        
        % difference for n1 
        diff_n1(t,b) = (L_current(b,1) - L_next(b,1))/(n1^0.9);
        
        % difference for n2
        diff_n2(t,b) = (L_current(b,2) - L_next(b,2))/(n2^0.9);
                     
    end
    L_current = L_next;
end

% looks for the smallest context tree such that the null hypothesis is rejected
pvalue = 1;
t = nTrees;
while (pvalue > alpha)&&(t > 1)
    t = t - 1;
    [~, pvalue] = ttest2(diff_n1(t,:), diff_n2(t,:), 'Alpha', alpha, 'Tail', 'right');
end
% if the null hypothesis was never rejected return the greatest tree
if pvalue > alpha, idx = 1; else, idx = t+1; end
% idx = t+1;
opt_tree = championTrees{idx};

end