function [contexts, Q, results] = estimate_discreteSeqROCTM(X, Y, Alphabet, varargin)
%ESTIMATE_DISCRETESEQROCTM Model selection for a SeqROCTM (X,Y)
%   [CONTEXTS, Q] = ESTIMATE_DISCRETESEQROCTM(X,Y,A) estimates a context
%   tree and a family of distributions Q (parameters of a SeqROCTM model)
%   from the sequences X,Y taking values in the alphabet A. The estimate
%   context tree is returned in CONTEXT and the family of distributions is
%   returned in Q.
%
%   [CONTEXTS, Q, RESULTS] = ESTIMATE_DISCRETESEQROCTM(...) returns a structure
%   with the following fields (only when 'bic' is chosen):
%       'nodes'             --  all possible nodes of the tree 
%       'stats'             --  statistics of the nodes [log_Ps, \sum log_Vas, log_Vs, \Xi_s, Nw, Nwa]
%       'nonExistingNodes'  --  nodes of the tree that do not appear in the
%                               sequence X
%       'XlengthWithoutNaN' --  length of the sequence X (because this
%                               function does not treat Nan values)
%       'nonanIndexes'      --  indexes of the elements in the sequence X
%                               (because this function does not treat Nan
%                               values)
%
%   [...] = ESTIMATE_DISCRETESEQROCTM(X,A,'PARAM1',val1,'PARAM2',val2,...)
%   specifies one or more of the following name/value pairs:
%
%       Parameter                Value
%       'EstimationMethod'       'bic' to estimate the model using the
%                                Bayesian Information Criteria.
%                                'context_cL' to estimate the model using the
%                                Context Algorithm based on likelihoods.
%                                'context_empD' to estimate model using the
%                                algorithm Context based on the comparison
%                                of distributions.
%                                Default is 'bic'.
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
%                                position the length of X (used to speed-up).
%                                Default value is [].

%   References:
%      [1] N. Hernández et al., arXiv xxx, (2021).   

%
%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 01/2021


%%%%%%%% name-value pairs arguments
% default values
options = struct('EstimationMethod', 'bic', 'MaxTreeHeight', floor(log(length(X))), ...
                    'ParameterValue', 1, 'CtxCompleteTree', -1, 'CtxTestStructure', -1, 'BicDegreeOfFreedom', 'fix', ...
                        'BicPrecomputedStats', []);

% acceptable names
optionNames = fieldnames(options);

for pair = reshape(varargin, 2, [])
    inpName = pair{1};
    
    if any(strcmpi(inpName, optionNames))
        options.(inpName) = pair{2};
    else
        error('%s is not a recognized parameter name', inpName);
    end
end
%%%%%%%%%%%%%%%%%%

if strcmpi('bic', options.EstimationMethod)
    df = ~strcmpi(options.BicDegreeOfFreedom, 'fix');
    if isempty(options.BicPrecomputedStats)
        [contexts, Q, results] = bic_WCT(X, Alphabet, options.MaxTreeHeight, options.ParameterValue, df, 0, Y);
    else
        [contexts, Q, results] = bic_WCT(X, Alphabet, options.MaxTreeHeight, options.ParameterValue, df, 0, Y, options.BicPrecomputedStats);
    end
elseif any(strcmpi(options.EstimationMethod, {'context_empD','context_cL'}))
    if isequal(options.CtxCompleteTree, -1)
        [contexts, Q] = CTestimator(X, Alphabet, options.MaxTreeHeight, options.EstimationMethod, options.ParameterValue, Y);
    else
        if isequal(options.CtxTestStructure, -1)
            [contexts, Q] = CTestimator(X, Alphabet, options.MaxTreeHeight, options.EstimationMethod, options.ParameterValue, Y, options.CtxCompleteTree);
        else
            [contexts, Q] = CTestimator(X, Alphabet, options.MaxTreeHeight, options.EstimationMethod, options.ParameterValue, Y, options.CtxCompleteTree, options.CtxTestStructure);
        end
    end
else
    error('%s is not a recognized parameter value', options.EstimationMethod);
end

end
