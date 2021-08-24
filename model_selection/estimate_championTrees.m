function [Trees, P, ML, cutoff] = estimate_championTrees(X, A, varargin)
%ESTIMATE_CHAMPIONTREES Compute the Champion Trees.
%   [TREE, P] = ESTIMATE_CHAMPIONTREES(X, A) compute the champion trees
%   using the sequence X taking values in the alphabet A. The champion
%   trees are returned in the cell array TREE and the corresponding family
%   of distributions in the cell array P.
%
%   [TREE, P, ML, CUTOFF] = ESTIMATE_CHAMPIONTREES(...) returns in the
%   vector ML the likelihood of each champion tree and in vector CUTOFF the
%   value of the hyperparameter with which it was obtained each champion
%   tree.
%
%   [...] = ESTIMATE_CHAMPIONTREES(X,A,'PARAM1',val1,'PARAM2',val2,...)
%   specifies one or more of the following name/value pairs:
%
%       Parameter                Value
%       'EstimationMethod'       'bic' to estimate the context tree models
%                                using the Bayesian Information Criteria or
%                                'context_cL' to estimate the context tree
%                                models using the Context Algorithm base on
%                                comparison of likelihoods.
%                                'context_empD' to estimate the context tree
%                                models using the Context Algorithm based
%                                on comparison of distributions.
%                                Default is 'bic'.
%       'MaxTreeHeight'          Maximum height of the context tree.
%                                Default is log(length(X)).
%       'ParameterLowerBound'    Minimum value of the parameter to be
%                                tuned. Default is 0.
%       'ParameterUpperBound'    Maximum value of the parameter to be
%                                tuned. Default is 100.
%       'Tolerance'              Minimum distance between parameter values.
%                                Default value is 10^-5.
%       'BicDegreeOfFreedom'     Degree of freedom used during the
%                                penalization in the BIC algorithm. 'fix'
%                                => (|A|-1), 'variable' => \sum_{a \in A}
%                                1{P(a|w)~=0}. Default value is 'fix'.
%       'BicMissing'             0 if there are no missing values in the
%                                sample, 1 is there are missing values.
%                                Default value is 0.

%   References:
%      [1] A. Galves et al., Ann. Appl. Stat., 6, 1, 186-209 (2012)
%      [2] N. Hernández et al., arXiv:2009.06371, (2021).   

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 01/2021

%%%%%%%% name-value pairs arguments
% default values
options = struct('EstimationMethod', 'bic', 'MaxTreeHeight', floor(log(length(X))), ...
                    'ParameterLowerBound', 0, 'ParameterUpperBound', 100, 'Tolerance', 10^-5,...
                    'BicDegreeOfFreedom', 'fix', 'BicMissing', 0);

% acceptable names
optionNames = fieldnames(options);

for pair = reshape(varargin, 2, [])
    inpName = pair{1};
    
    if any(strcmp(inpName, optionNames))
        
        if strcmp(inpName, 'EstimationMethod')
            if any(strcmpi(pair{2}, {'bic','context_empD', 'context_cL'}))
                options.(inpName) = pair{2};
            else
                error('%s is not a recognized parameter value', pair{2})
            end
        else
            options.(inpName) = pair{2};
        end
        
    else
        error('%s is not a recognized parameter name', inpName);
    end
end
%%%%%%%%%%%%%%%%%%

l_min = options.ParameterLowerBound;
u = options.ParameterUpperBound;
tol = options.Tolerance;
max_height = options.MaxTreeHeight;
df = options.BicDegreeOfFreedom;
miss = options.BicMissing;

lA = length(A);
if any(strcmpi(options.EstimationMethod, {'context_empD','context_cL'}))
    % compute the complete tree and the TEST structure only once (for speed-up)
    [T, I] = completetree(X, max_height, A);
    TEST = getTESTstructure(T, I, lA, X);
    ct_inf = {T, I};
    
    %BIC-info
    precomputed_stats = [];
else
    % compute the statistics Nw and Nwa used in BIC only once (for speed-up)
    df1 = ~strcmpi(options.BicDegreeOfFreedom,'fix');
    [~, ~, ~, outps] = bic_WCT(X, A, max_height, l_min, df1, options.BicMissing);
    precomputed_stats{1} = outps.stats(:,5:6+lA-1);
    precomputed_stats{2} = outps.nonExistingNodes;
    precomputed_stats{3} = outps.XlengthWithoutNaN;
    nonanIndexes = outps.nonanIndexes;
    
    %ctx-info
    TEST = -1;
    ct_inf = -1;
end

% if treatment of missing values is active, compute the indexes that are
% going to be passed to the likelihood function
if options.BicMissing
    param_likhd = nonanIndexes;
else
    param_likhd = 0;
end

% estimate the trees for the minimum and maximal value of the penalization
% constant
[tau_l, p_l] = estimate_contexttree(X, A, 'MaxTreeHeight', max_height, 'EstimationMethod', options.EstimationMethod, 'ParameterValue', l_min, 'CtxCompleteTree', ct_inf, 'CtxTestStructure', TEST, 'BicDegreeOfFreedom', df, 'BicPrecomputedStats', precomputed_stats, 'BicMissing', miss);  
[tau_upper, p_upper] = estimate_contexttree(X, A, 'MaxTreeHeight', max_height, 'EstimationMethod', options.EstimationMethod, 'ParameterValue', u, 'CtxCompleteTree', ct_inf, 'CtxTestStructure', TEST, 'BicDegreeOfFreedom', df, 'BicPrecomputedStats', precomputed_stats, 'BicMissing', miss);

% initialize
upper_bound = u;
Trees = {}; P = {}; ML = []; cutoff = [];

if ~isempty(tau_upper)
    disp('Warning: The empty tree is not obtain for the maximum value of the penalization constant given.')
end

    tau_u = tau_upper;
    p_u = p_upper;

    i = 1; Trees{i} = tau_l; P{i} = p_l; ML(i) = treeloglikelihood(X, tau_l, A, param_likhd); cutoff(i) = l_min;

    % estimate the different trees in the interval specified for the
    % penalization constant
    while ~isequalCT(tau_l, tau_upper)
        while abs(u - l_min) > tol
            while ~isequalCT(tau_u, tau_l)&&(abs(u - l_min) > tol) % the second condition its necessary because for some cases,
                a = u;                                               % the complete tree is obtain when l_min=0 and for any value     
                tau_a = tau_u; p_a = p_u;                            % greater than zero, a tree different from the complete tree is obtained   
                u = (l_min + u)/2;                                      
                [tau_u, p_u] = estimate_contexttree(X, A, 'MaxTreeHeight', max_height, 'EstimationMethod', options.EstimationMethod,...
                                                    'ParameterValue', u, 'CtxCompleteTree', ct_inf, 'CtxTestStructure', TEST, ...
                                                    'BicDegreeOfFreedom', df, 'BicPrecomputedStats', precomputed_stats, ...
                                                    'BicMissing', miss);   
            end
            l_min = u; tau_l = tau_u;
            u = a; tau_u = tau_a; p_u = p_a;
        end
        i = i + 1;
        Trees{i} = tau_u; P{i} = p_u; cutoff(i) = u;
        ML(i) = treeloglikelihood(X, tau_u, A, param_likhd);
        l_min = u; tau_l = tau_u; 
        u = upper_bound;
        tau_u = tau_upper;  
        p_u = p_upper;
    end
end