function [opt_tree, idx, bootsamples_n2] = modeltunning_SMC(championTrees, A, n1, n2, alpha, B, bootstrategy, param1, param2)
%MODELTUNNING_SMC Context tree selection using the Smallest Maximizer Criterion 
%				  (see article Galves, A. et. al., Ann. Appl. Stat., Volume 6, Number 1 (2012), 186-209)
%
% Inputs
%
%   championTrees : set of Champion Trees (trees obtained for different
%                   values of the penalization constant in the BIC
%                   criteria)
%   A             : alphabet
%   n1            : proportion of the size of the sample corresponding to
%                   the size of the smaller resample.
%   n2            : proportion of the size of the sample corresponding to
%                   the size of the larger resample.
%   alpha         : alpha level to use on the t-test
%   B             : number of resamples in the bootstrap procedure
%   bootstrategy  : bootstrap procedure used. 'parametric_ctm': the largest
%                   tree in the set of Champion Trees is used to generate
%                   the bootstrap samples. 'blocks': a renewal point is used
%                   to create blocks and use them to generate the bootstrap
%                   samples.
%  param1, param2 : when the bootstrap strategy is 'parametric_ctm' refer to
%                   the context tree model (tree and distributions) used to
%                   generate the samples. When the bootstrap strategy is
%                   'blocks' refer to the sequence and the renewal point.
%
% Outputs
%
%  opt_tree         : optimal context tree
%  idx              : index of the optimal context tree in the set of
%                       Champion Trees
%  bootsamples_n2   : bootstrap samples generated
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 07/2020

nTrees = length(championTrees);

% generate B bootstrap samples of size n2
switch bootstrategy
    case 'blocks'
        X = param1;
        renewal_point = param2;
        bootsamples_n2 = bootstrap_blocks(X, renewal_point, n2, B);
    case 'parametric_ctm'
        tau0 = param1;
        P0 = param2;
        bootsamples_n2 = generatesampleCTM_fast(tau0, P0, A, n2, B);
end

% compute the differences in likelihood for each pair of consecutive trees
% and all the bootstrap samples
diff_n1 = zeros(nTrees-1, B);
diff_n2 = zeros(nTrees-1, B);

% initialize L_current
L_current = zeros(B,2);
for b = 1 : B
    L_current(b,1) = treeloglikelihood(championTrees{1}, A, bootsamples_n2(b, 1:n1));
    L_current(b,2) = treeloglikelihood(championTrees{1}, A, bootsamples_n2(b,:));
end

for t = 1 : nTrees-1
    L_next = zeros(B,2); % store the log-likelihood of tree t+1 to speed-up
    for b = 1 : B
        %
        L_next(b,1) = treeloglikelihood(championTrees{t+1}, A, bootsamples_n2(b, 1:n1));
        L_next(b,2) = treeloglikelihood(championTrees{t+1}, A, bootsamples_n2(b,:));
        
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
idx = t+1;
opt_tree = championTrees{idx};

end