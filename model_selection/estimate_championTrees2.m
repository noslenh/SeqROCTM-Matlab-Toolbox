function [Trees, P, ML, cutoff] = estimate_championTrees2(X, Y, A, varargin)
%ESTIMATE_CHAMPIONTREES estimate the Champion Trees for the SeqROCTM (X,Y).
%                       Using the 'bic' criterion a context tree is
%                       estimated for different values of the penalization
%                       constant. This function implements the method
%                       proposed in Ann. Appl. Stat. 6 (1), 2012, 186-209
%                       for SeqROCTM
%
% Input
%   
%   X           : sequence of inputs
%   Y           : sequence of responses
%   max_height  : height of the complete tree
%   A           : alphabet
%   varargin    : contain the l_min and u values
%                 l_min: minimum value for the penalization constant (usually 0)
%                 u    : maximum value for the penalization constant
%                        (usually a big value such that the resulting tree is the empty tree)
%                 tol  : tolerance used when estimating the set of champion trees  
%
% Output
%   
%   Trees       : set of champion trees
%   P           : set of distributions associated to the trees
%   ML          : maximum likelihood value for each of the champion tree
%   cuttoff     : values of the bic penalization
%
%
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 01/2021

%%%%%%%% name-value pairs arguments
% default values
options = struct('EstimationMethod', 'bic', 'MaxTreeHeight', log(length(X)), ...
    'ParameterLowerBound', 0, 'ParameterUpperBound', 100, 'Tolerance', 10^-5, 'DegreeOfFreedom', 'fix');

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
df = options.DegreeOfFreedom;

lA = length(A);
lX = length(X);

if strcmpi(options.EstimationMethod, 'context')
    % compute the complete tree and the TEST structure only once (for speed-up)
    [T, I] = completetree(X, max_height, A);
    TEST = getTESTstructure(T, I, lA, Y);
    ct_inf = {T, I};
else
    %
    TEST = -1;
    ct_inf = -1;
    % compute the statistics Nw and Nwa used in BIC only once (for speed-up)
    df1 = ~strcmpi(options.DegreeOfFreedom,'fix');
    [~,~,~,~,~,~,~, stats, non_existing_nodes] = get_maximizingTree([], lA, max_height, [], X, lX, l_min, df1, 0, Y);
    precomputed_stats{1} = stats(:,5:6+lA-1);
    precomputed_stats{2} = non_existing_nodes;
end

% estimate the trees for the minimum and maximal value of the penalization
% constant
[tau_l, p_l] = estimate_discreteSeqROCTM(X, Y, A, 'MaxTreeHeight', max_height, 'EstimationMethod', options.EstimationMethod, 'ParameterValue', l_min, 'CompleteTree', ct_inf, 'TestStructure', TEST, 'DegreeOfFreedom', df, 'BicPrecomputedStats', precomputed_stats);  
[tau_upper, p_upper] = estimate_discreteSeqROCTM(X, Y, A, 'MaxTreeHeight', max_height, 'EstimationMethod', options.EstimationMethod, 'ParameterValue', u, 'CompleteTree', ct_inf, 'TestStructure', TEST, 'DegreeOfFreedom', df, 'BicPrecomputedStats', precomputed_stats);

% [tau_l, p_l] = CTestimator(X, A, max_height, 'context', l_min, Y, {T, I}, TEST);
% [tau_upper, p_upper] = CTestimator(X, A, max_height, 'context', u, Y, {T, I}, TEST);

% initialize
upper_bound = u;
Trees = {}; P = {}; ML = []; cutoff = [];

if ~isempty(tau_upper)
    disp('Warning: The empty tree is not obtain for the maximun value of the penalization constant given.')
end

    tau_u = tau_upper;
    p_u = p_upper;

    i = 1; Trees{i} = tau_l; P{i} = p_l; ML(i) = treeloglikelihood(tau_l, A, X, Y); cutoff(i) = l_min;

    % estimate the different trees in the interval specified for the
    % penalization constant
    while ~isequalCT(tau_l, tau_upper)
        while abs(u - l_min) > tol
            while ~isequalCT(tau_u, tau_l)&&(abs(u - l_min) > 10^-5) % the second condition its necessary because for some cases,
                a = u;                                               % the complete tree is obtain when l_min=0 and for any value     
                tau_a = tau_u; p_a = p_u;                            % greater than zero, a tree different from the complete tree is obtained   
                u = (l_min + u)/2;                                      
%                 [tau_u, p_u] = CTestimator(X, A, max_height, 'context', u, Y, {T, I}, TEST);  
                [tau_u, p_u] = estimate_discreteSeqROCTM(X, Y, A, 'MaxTreeHeight', max_height, 'EstimationMethod', options.EstimationMethod, 'ParameterValue', u, 'CompleteTree', ct_inf, 'TestStructure', TEST, 'DegreeOfFreedom', df, 'BicPrecomputedStats', precomputed_stats);
            end
            l_min = u; tau_l = tau_u;
            u = a; tau_u = tau_a; p_u = p_a;
        end
        i = i + 1;
        Trees{i} = tau_u; P{i} = p_u; cutoff(i) = u;
        ML(i) = treeloglikelihood(tau_u, A, X, Y);
        l_min = u; tau_l = tau_u; 
        u = upper_bound;
        tau_u = tau_upper;  
        p_u = p_upper;
    end
end