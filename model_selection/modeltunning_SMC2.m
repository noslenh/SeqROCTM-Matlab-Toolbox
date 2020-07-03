function [opt_tree, idx, Xbootsamples_n2, Ybootsamples_n2] = modeltunning_SMC2(championTrees, A, n1, n2, alpha, B, ... 
                                                                bootstrategyX, param1, param2, tau_y, q_y)
%MODELTUNNING_SMC Context tree selection using the Smallest Maximizer
%                 Criterion for SeqROCTM 
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
%   Xbootstrategy : bootstrap procedure used. 'parametric_ctm': a given
%                   context tree model is used to generate the boostrap
%                   samples. 'blocks': a renewal point is used to create
%                   independent blocks that are used to generate the bootstrap
%                   samples.
%  param1, param2 : when the bootstrap strategy is 'parametric_ctm' refer to
%                   the context tree model (tree and distributions) used to
%                   generate the samples. When the bootstrap strategy is
%                   'blocks' refer to the sequence and the renewal point.
%  tau_y           : context tree used to generate the bootstrap samples for
%                   the response sequence.
%  q_y             : distributions associated to the contexts in tau_y  
%
% Outputs
%
%  opt_tree         : optimal context tree
%  idx              : index of the optimal context tree in the set of
%                       Champion Trees
%  Xbootsamples_n2  : bootstrap samples used for the X sequence
%  Ybootsamples_n2  : bootstrap samples used for the Y sequence

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 07/2020

nTrees = length(championTrees);
ml = max(cellfun(@(x) length(x), tau_y));

% generate B bootstrap samples of size n2 of the sequence X
switch bootstrategyX
    case 'blocks'
        X = param1;
        renewal_point = param2;
        Xboot = bootstrap_blocks(X, renewal_point, n2+ml, B);
    case 'parametric_ctm'
        tau0 = param1;
        P0 = param2;
        Xboot = generatesampleCTM_fast(tau0, P0, A, n2+ml, B);
    case 'none'
        X = param1;
        lX = size(X,2);
        if n2+ml <= lX
            Xboot = ones(B,1) * X(:, 1:n2+ml);
        else
            disp(['The n2 parameter must be lower than ' num2str(lX - ml)]);
            return
        end
end

% 'parametric_ctm' is the only way (in this version) to bootstrap the
% response sequence
% generate B bootstrap samples of size n2 of the sequence Y
        
Ybootsamples_n2 = zeros(B, n2);
Xbootsamples_n2 = zeros(B, n2);
for i = 1 : B
    [Xbootsamples_n2(i,:), Ybootsamples_n2(i,:)] = generatesampleYSeqROCTM(Xboot(i,:), tau_y, q_y, A);
end


% compute the differences in likelihood for each pair of consecutives trees
% and all the bootstrap samples
diff_n1 = zeros(nTrees-1, B);
diff_n2 = zeros(nTrees-1, B);

% initialize L_current
L_current = zeros(B,2);
for b = 1 : B
    L_current(b,1) = treeloglikelihood(championTrees{1}, A, Xbootsamples_n2(b, 1:n1), Ybootsamples_n2(b, 1:n1));
    L_current(b,2) = treeloglikelihood(championTrees{1}, A, Xbootsamples_n2(b,:), Ybootsamples_n2(b,:));
end

for t = 1 : nTrees-1
    L_next = zeros(B,2); % store the log-likelihood of tree t+1 to speed-up
    for b = 1 : B
        %
        L_next(b,1) = treeloglikelihood(championTrees{t+1}, A, Xbootsamples_n2(b, 1:n1), Ybootsamples_n2(b, 1:n1));
        L_next(b,2) = treeloglikelihood(championTrees{t+1}, A, Xbootsamples_n2(b,:), Ybootsamples_n2(b,:));
        
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