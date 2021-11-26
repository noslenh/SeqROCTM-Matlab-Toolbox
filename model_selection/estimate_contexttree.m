function [contexts, P, results] = estimate_contexttree(X, A, varargin)
%ESTIMATE_CONTEXTTREE Estimate a context tree from the sequence X.
%   [CONTEXTS, P] = ESTIMATE_CONTEXTTREE(X,A) estimates a context tree
%   model from the sequence X taking values in the alphabet A. The
%   estimated context tree is returned in CONTEXTS and the corresponding
%   family of distributions is returned in P.
%
%   [CONTEXTS, P, RESULTS] = ESTIMATE_CONTEXTTREE(...) returns a structure
%   with the following fields (only when 'bic' is chosen):
%       'nodes'             --  all possible nodes of the tree 
%       'stats'             --  statistics of the nodes [log_Ps, \sum log_Vas, log_Vs, \Xi_s, Nw, Nwa]
%       'nonExistingNodes'  --  nodes of the tree that do not appear in the
%                               sequence X
%       'XlengthWithoutNaN' --  length of the sequence X discarding the Nan
%                               values
%       'nonanIndexes'      --  indexes of the non Nan elements in the
%                               sequence X
%
%   [...] = ESTIMATE_CONTEXTTREE(X,A,'PARAM1',val1,'PARAM2',val2,...)
%   specifies one or more of the following name/value pairs:
%
%       Parameter                Value
%       'EstimationMethod'       'bic' to estimate the context tree model
%                                using the Bayesian Information Criteria.
%                                'context' to estimate the context tree
%                                model using the Context Algorithm.
%                                'emp_distribution' to estimate the context
%                                tree model using a pruning criterion based
%                                on the comparison of the empirical
%                                transition probabilities associated to the
%                                leaves of the tree. Default is 'bic'.
%       'MaxTreeHeight'          Maximum height of the context tree.
%                                Default is log(length(X)).
%       'ParameterValue'         Value of the hyperparameter involved on
%                                each estimation method: the penalization
%                                constant in the case of 'bic' and the
%                                threshold in the case of 'context' and
%                                'emp_distribution'. Default value is 1.
%       'CtxCompleteTree'        Initial complete tree for the algorithm
%                                Context (usually used to speed-up).
%                                Default value is [].
%       'CtxTestStructure'       Auxiliary structure for the algorithm
%                                Context (usually used to speed-up).
%                                Default value is [].
%       'BicDegreeOfFreedom'     Degree of freedom used during the
%                                penalization in the BIC algorithm. 'fix'
%                                => (|A|-1), 'variable' => \sum_{a \in A}
%                                1{Q(a|w)~=0}. Default value is 'fix'.
%       'BicPrecomputedStats'    A cell array containing in the first
%                                position [Nw Nwa] for each node. In the
%                                second position the nodes that do not
%                                appear in the sequence and in the last
%                                position the length of X without NaN
%                                (usually used to speed-up). Default value
%                                is [].
%       'BicMissing'             1 if it is required treatment of Nan
%                                values, 0 otherwise

%   References:
%      [1] J. Rissanen, IEEE Trans. Inform. Theory 29, 656-664 (1983)
%      [2] I. Csiszar et al., IEEE Trans. Inform. Theory, 3, 52, 1007-1016 (2006)
%      [3] A. Galves et al., Progress in Probability 60, 257-269 (2008)

%
%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 01/2021

%%%%%%%% name-value pairs arguments
% default values
options = struct('EstimationMethod', 'bic', 'MaxTreeHeight', floor(log(length(X))), ...
                    'ParameterValue', 1, 'CtxCompleteTree', -1, 'CtxTestStructure', -1, 'BicDegreeOfFreedom', 'fix', ...
                        'BicPrecomputedStats', [], 'BicMissing', 0);

% acceptable names
optionNames = fieldnames(options);

for pair = reshape(varargin, 2, [])
    inpName = pair{1};
    
    if any(strcmp(inpName, optionNames))
        options.(inpName) = pair{2};
    else
        error('%s is not a recognized parameter name', inpName);
    end
end
%%%%%%%%%%%%%%%%%%

if strcmpi('bic', options.EstimationMethod)
    df = ~strcmpi(options.BicDegreeOfFreedom,'fix');
    if isempty(options.BicPrecomputedStats)
        [contexts, P, ~, results] = bic_WCT(X, A, options.MaxTreeHeight, options.ParameterValue, df, options.BicMissing);
    else
        [contexts, P, ~, results] = bic_WCT(X, A, options.MaxTreeHeight, options.ParameterValue, df, options.BicMissing, [], options.BicPrecomputedStats);
    end
elseif any(strcmpi(options.EstimationMethod, {'context','emp_distribution'}))
    if isequal(options.CtxCompleteTree, -1)
        [contexts, P] = CTestimator(X, A, options.MaxTreeHeight, options.EstimationMethod, options.ParameterValue);
    else
        if isequal(options.CtxTestStructure, -1)
            [contexts, P] = CTestimator(X, A, options.MaxTreeHeight, options.EstimationMethod, options.ParameterValue, [], options.CtxCompleteTree);
        else
            [contexts, P] = CTestimator(X, A, options.MaxTreeHeight, options.EstimationMethod, options.ParameterValue, [], options.CtxCompleteTree, options.CtxTestStructure);
        end
    end
else
    error('%s is not a recognized parameter value', options.EstimationMethod);
end

end
