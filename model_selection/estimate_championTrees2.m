function [Trees, Q, ML, cutoff] = estimate_championTrees2(X, Y, A, varargin)
%ESTIMATE_CHAMPIONTREES Compute the Champion Trees.
%   [TREE, P] = ESTIMATE_CHAMPIONTREES(X, A) compute the champion trees
%   using the SeqROCTM (X,Y) taking values in the alphabet A. The champion
%   trees are returned in the cell array TREE and the corresponding family
%   of distributions in the cell array Q.
%
%   [TREE, Q, ML, CUTOFF] = ESTIMATE_CHAMPIONTREES(...) returns in the
%   vector ML the likelihood of each champion tree and in vector CUTOFF the
%   value of the hyperparameter with wich it was obtained each champion
%   tree.
%
%   [...] = ESTIMATE_CHAMPIONTREES(X,A,'PARAM1',val1,'PARAM2',val2,...)
%   specifies one or more of the following name/value pairs:
%
%       Parameter                Value
%       'EstimationMethod'       'bic' to estimate the context tree models
%                                using the Bayesian Information Criteria or
%                                'context' to estimate the context tree
%                                models using the Context Algorithm.
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
%                                1{Q(a|w)~=0}. Default value is 'fix'.

%   References:
%      [1] A. Galves et al., Ann. Appl. Stat., 6, 1, 186-209 (2012)
%      [2] N. Hernández et al., arXiv xxx, (2021).   

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 01/2021


%%%%%%%% name-value pairs arguments
% default values
options = struct('EstimationMethod', 'bic', 'MaxTreeHeight', log(length(X)), ...
    'ParameterLowerBound', 0, 'ParameterUpperBound', 100, 'Tolerance', 10^-5, 'BicDegreeOfFreedom', 'fix');

% acceptable names
optionNames = fieldnames(options);

for pair = reshape(varargin, 2, [])
    inpName = pair{1};
    
    if any(strcmpi(inpName, optionNames))
        
        if strcmpi(inpName, 'estimationmethod')
            if any(strcmpi(pair{2}, {'bic','context'}))
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

lA = length(A);
if strcmpi(options.EstimationMethod, 'context')
    % compute the complete tree and the TEST structure only once (for speed-up)
    [T, I] = completetree(X, max_height, A);
    TEST = getTESTstructure(T, I, lA, Y);
    ct_inf = {T, I};
    
    %BIC-info
    precomputed_stats = [];
else
    % compute the statistics Nw and Nwa used in BIC only once (for speed-up)
    df1 = ~strcmpi(options.BicDegreeOfFreedom,'fix');
    [~, ~, ~, outps] = bic_WCT(X, A, max_height, l_min, df1, 0, Y);
    precomputed_stats{1} = outps.stats(:,5:6+lA-1);
    precomputed_stats{2} = outps.nonExistingNodes;
    precomputed_stats{3} = outps.XlengthWithoutNaN;
    
    % ctx-info
    TEST = -1;
    ct_inf = -1;
end

% estimate the trees for the minimum and maximal value of the penalization
% constant
[tau_l, p_l] = estimate_discreteSeqROCTM(X, Y, A, 'MaxTreeHeight', max_height, 'EstimationMethod', options.EstimationMethod, 'ParameterValue', l_min, 'CtxCompleteTree', ct_inf, 'CtxTestStructure', TEST, 'BicDegreeOfFreedom', df, 'BicPrecomputedStats', precomputed_stats);  
[tau_upper, p_upper] = estimate_discreteSeqROCTM(X, Y, A, 'MaxTreeHeight', max_height, 'EstimationMethod', options.EstimationMethod, 'ParameterValue', u, 'CtxCompleteTree', ct_inf, 'CtxTestStructure', TEST, 'BicDegreeOfFreedom', df, 'BicPrecomputedStats', precomputed_stats);

% [tau_l, p_l] = CTestimator(X, A, max_height, 'context', l_min, Y, {T, I}, TEST);
% [tau_upper, p_upper] = CTestimator(X, A, max_height, 'context', u, Y, {T, I}, TEST);

% initialize
upper_bound = u;
Trees = {}; Q = {}; ML = []; cutoff = [];

if ~isempty(tau_upper)
    disp('Warning: The empty tree is not obtain for the maximun value of the penalization constant given.')
end

    tau_u = tau_upper;
    p_u = p_upper;

    i = 1; Trees{i} = tau_l; Q{i} = p_l; ML(i) = treeloglikelihood2(X, Y, tau_l, A); cutoff(i) = l_min;

    % estimate the different trees in the interval specified for the
    % penalization constant
    while ~isequalCT(tau_l, tau_upper)
        while abs(u - l_min) > tol
            while ~isequalCT(tau_u, tau_l)&&(abs(u - l_min) > tol) % the second condition its necessary because for some cases,
                a = u;                                               % the complete tree is obtain when l_min=0 and for any value     
                tau_a = tau_u; p_a = p_u;                            % greater than zero, a tree different from the complete tree is obtained   
                u = (l_min + u)/2;                                      
%                 [tau_u, p_u] = CTestimator(X, A, max_height, 'context', u, Y, {T, I}, TEST);  
                [tau_u, p_u] = estimate_discreteSeqROCTM(X, Y, A, 'MaxTreeHeight', max_height, 'EstimationMethod', options.EstimationMethod, 'ParameterValue', u, 'CtxCompleteTree', ct_inf, 'CtxTestStructure', TEST, 'BicDegreeOfFreedom', df, 'BicPrecomputedStats', precomputed_stats);
            end
            l_min = u; tau_l = tau_u;
            u = a; tau_u = tau_a; p_u = p_a;
        end
        i = i + 1;
        Trees{i} = tau_u; Q{i} = p_u; cutoff(i) = u;
        ML(i) = treeloglikelihood2(X, Y, tau_u, A);
        l_min = u; tau_l = tau_u; 
        u = upper_bound;
        tau_u = tau_upper;  
        p_u = p_upper;
    end
end