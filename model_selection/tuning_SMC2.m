function [opt_tree, idx] = tuning_SMC2(championTrees, A, n1, n2, alpha, Xbootsamples, Ybootsamples)
%MODELTUNNING_SMC Context tree selection using the Smallest Maximizer
%                 Criterion for SeqROCTM 
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
%   alpha         : significance level used in the t-test
%   Xbootsamples  : bootstrap samples for the X sequence
%   Ybootsamples  : bootstrap samples for the Y sequence
%
% Outputs
%
%  opt_tree         : optimal context tree
%  idx              : index of the optimal context tree in the set of
%                       Champion Trees

%   References:
%      [1] A. Galves et al., Ann. Appl. Stat., 6, 1, 186-209 (2012)
%      [2] N. Hernández et al., arXiv xxx, (2021).   

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 01/2021

nTrees = length(championTrees);
B = size(Xbootsamples,1);

% ml = max(cellfun(@(x) length(x), tau_y));

% % generate B bootstrap samples of size n2 of the sequence X
% switch bootstrategyX
%     case 'blocks'
%         X = param1;
%         renewal_point = param2;
%         Xboot = bootstrap_blocks(X, renewal_point, n2+ml, B);
%     case 'parametric_ctm'
%         tau0 = param1;
%         P0 = param2;
%         Xboot = generatesampleCTM_fast(tau0, P0, A, n2+ml, B);
%     case 'none'
%         X = param1;
%         lX = size(X,2);
%         if n2+ml <= lX
%             Xboot = ones(B,1) * X(:, 1:n2+ml);
%         else
%             disp(['The n2 parameter must be lower than ' num2str(lX - ml)]);
%             return
%         end
% end
% 
% % 'parametric_ctm' is the only way (in this version) to bootstrap the
% % response sequence
% % generate B bootstrap samples of size n2 of the sequence Y
%         
% Ybootsamples_n2 = zeros(B, n2);
% Xbootsamples_n2 = zeros(B, n2);
% for i = 1 : B
%     [Xbootsamples_n2(i,:), Ybootsamples_n2(i,:)] = generatesampleYSeqROCTM(Xboot(i,:), tau_y, q_y, A);
% end
% 

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
% idx = t+1;
opt_tree = championTrees{idx};

end