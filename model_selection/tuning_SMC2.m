function [opt_tree, idx] = tuning_SMC2(championTrees, A, n1, n2, alpha, Xbootsamples, Ybootsamples)
%MODELTUNNING_SMC Context tree selection using the Smallest Maximizer Criterion for SeqROCTM.
%
% Inputs
%
%   championTrees : set of Champion Trees (trees obtained for different
%                   values of the penalization constant in the BIC
%                   criteria)
%   A             : alphabet
%   n1            : proportion of the size of the sample corresponding to
%                   the size of the smaller re-sample.
%   n2            : proportion of the size of the sample corresponding to
%                   the size of the larger re-sample.
%   alpha         : significance level used in the t-test
%   Xbootsamples  : bootstrap samples for the X sequence
%   Ybootsamples  : bootstrap samples for the Y sequence
%
% Outputs
%
%  opt_tree         : optimal context tree
%  idx              : index of the optimal context tree in the set of
%                       Champion Trees
%
%   References:
%      [1] A. Galves et al., Ann. Appl. Stat., 6, 1, 186-209 (2012)
%      [2] N. Hernández et al., arXiv:2009.06371, (2021).  
%
%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 01/2021

nTrees = length(championTrees);
B = size(Xbootsamples,1);

% compute the differences in likelihood for each pair of consecutive trees
% and all the bootstrap samples
diff_n1 = zeros(nTrees-1, B);
diff_n2 = zeros(nTrees-1, B);

% initialize L_current
L_current = zeros(B,2);
for b = 1 : B
    L_current(b,1) = treeloglikelihood2(Xbootsamples(b, 1:n1), Ybootsamples(b, 1:n1), championTrees{1}, A);
    L_current(b,2) = treeloglikelihood2(Xbootsamples(b,:), Ybootsamples(b,:), championTrees{1}, A);
end

for t = 1 : nTrees-1
    L_next = zeros(B,2); % store the log-likelihood of tree t+1 to speed-up
    for b = 1 : B
        %
        L_next(b,1) = treeloglikelihood2(Xbootsamples(b, 1:n1), Ybootsamples(b, 1:n1), championTrees{t+1}, A);
        L_next(b,2) = treeloglikelihood2(Xbootsamples(b,:), Ybootsamples(b,:), championTrees{t+1}, A);
        
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
opt_tree = championTrees{idx};

end