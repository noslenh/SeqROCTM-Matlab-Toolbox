function [optTree, optP, results] = tune_contextTreeModel(X, A, varargin)
%TUNE_CONTEXTTREEMODEL Tune a context tree estimation algorithm.
%   [OPTTREE, OPTP] = TUNE_CONTEXTTREEMODEL(X,A) tunes a context tree model
%   using the sequence X of values in the alphabet A. The optimal context
%   tree is returned in OPTTREE and the probability distributions
%   associated to the contexts are returned in OPTP.
%
%   [OPTTREE, OPTP, RESULTS] = TUNE_CONTEXTTREEMODEL(...) returns a structure with the following fields:
%       'champions'     -- the Champion Trees
%       'Ps'            -- the family of probability distributions for each champion tree
%       'fvalues'       -- likelihood or risk function values for the champion trees
%       'prmvalues'     -- the parameter value 
%       'idxOptTree'    -- index of the optimal model in the champion trees
%       'bootsamples'   -- bootstrap samples generated to tune the model
%
%   TUNE_CONTEXTTREEMODEL treats NaNs as missing values. The computational
%   cost of the algorithm can increase. 
%
%   [...] = TUNE_CONTEXTTREEMODEL(X,A,'PARAM1',val1,'PARAM2',val2,...)
%   specifies one or more of the following name/value pairs:
%
%       Parameter                Value
%       'TuningMethod'          'smc' to perform the Smaller Maximizer
%                                Criteria or 'risk' to use a risk function.
%                                Default is 'smc'.
%       'EstimationMethod'       'bic' to estimate the context tree models
%                                using the Bayesian Information Criteria
%                                and tune the penalization constant or
%                                'context' to estimate the context tree
%                                models using the Context Algorithm and
%                                tune the threshold involved in the Context
%                                Algorithm. Default is 'bic'.
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
%       'n1'                     minimum sample length of bootstrap samples
%                                in SMC. Default value is floor(0.3*length(X)).
%       'n2'                     maximum sample length of bootstrap samples 
%                                in SMC. Default value is floor(0.9*length(X)).
%       'Alpha'                  Statistical significance of the t-test in
%                                SMC. Default value is 0.01.
%       'BootNSamples'           Number of bootstrap samples. Default is 200.
%       'BootStrategy'           bootstrap procedure. 'parametric': the
%                                largest tree in the set of Champion Trees
%                                is used to generate the bootstrap samples
%                                or a model given in 'BootModel'. 
%                                'blocks': a renewal point is used to
%                                create independent blocks and these blocks
%                                are sampled to create the bootstrap
%                                samples. Default value is 'blocks'.
%       'BootRenewalPoint'       the renewal point to be used with the
%                                bootstrap strategy 'blocks'. This can be a
%                                subsequence or the string 'compute'. When
%                                it equals 'compute' a renewal point is
%                                compute from the largest context tree in
%                                the Champion Trees. Default value is
%                                'compute'.
%       'BootModel'              context tree model used to generate the
%                                bootstrap samples when the bootstrap
%                                strategy is 'parametric'. A cell array
%                                that contains in BootModel{1} a context
%                                tree and in BootModel{2} the transition
%                                probabilities. Default value is [] (in
%                                this case the largest model in the
%                                champion trees is used).
%       

%   References:
%      [1] A. Galves et al., Ann. Appl. Stat., 6, 1, 186-209 (2012)
%      [2] N. Hernández et al., arXiv xxx, (2021).   
%      [3] P. Buhlmann et al., Ann. Inst. Statist. Math, 52, 1, 287-315 (2000)
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 01/2021

lX = length(X);

% name-value pairs arguments
% default values
options = struct(   'TuningMethod', 'smc',              ...
                    'EstimationMethod', 'bic',          ...  
                    'MaxTreeHeight', log(lX),           ...
                    'ParameterLowerBound', 0,           ...
                    'ParameterUpperBound', 100,         ...
                    'Tolerance', 10^-5,                 ...
                    'BicDegreeOfFreedom', 'fix',        ...
                    'BicMissing', 0,                    ...
                    'n1', floor(0.3*lX),                ...                        
                    'n2', floor(0.9*lX),                ...
                    'Alpha', 0.01,                      ...
                    'BootNSamples', 200,                ...
                    'BootStrategy', 'blocks',           ...
                    'BootRenewalPoint', 'compute',      ...
                    'BootModel', []                     ...
                    );

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

% compute the Champion Trees
[results.champions, results.Ps, results.fvalues, results.prmvalues] = estimate_championTrees(X, A, ...
                                                                'MaxTreeHeight', options.MaxTreeHeight,             ...
                                                                'ParameterLowerBound', options.ParameterLowerBound, ...
                                                                'ParameterUpperBound', options.ParameterUpperBound, ...
                                                                'EstimationMethod', options.EstimationMethod,       ...
                                                                'Tolerance', options.Tolerance,                     ...
                                                                'BicDegreeOfFreedom', options.BicDegreeOfFreedom,   ...
                                                                'BicMissing', options.BicMissing                    ...
                                                                );

smc = 1;
if strcmpi(options.TuningMethod, 'smc')
    lengthbootsamples = options.n2;
elseif strcmpi(options.TuningMethod, 'risk')
    smc = 0;
    lengthbootsamples = lX+1;
else
    error('%s is not a recognized tunning method', options.TuningMethod);
end
                                                                    
% generate the Bootstrap samples
switch options.BootStrategy
    case 'blocks'
        %compute a renewal point from the largest model in the champion
        %trees ('compute') or use the specified renewal point
        if strcmpi(options.BootRenewalPoint, 'compute')
            renewal_point = tree_renewalpoint(results.champions{1}, results.Ps{1}, A, X);
        else
            renewal_point = options.BootRenewalPoint;
        end
        results.bootsamples = bootstrap_blocks(X, renewal_point, lengthbootsamples, options.BootNSamples);
        bootstrap_missing = options.BicMissing;
    case 'parametric'
        %if no model was specified through 'BootModel' use the largest
        %model in the Champion Trees. Otherwise, use the model in 'BootModel'
        if isempty(options.BootModel)
            tau0 = results.champions{1};
            p0 = results.Ps{1};
        else
            tau0 = options.BootModel{1};
            p0 = options.BootModel{2};
        end
        results.bootsamples = generatesampleCTM_fast(tau0, p0, A, lengthbootsamples, options.BootNSamples);
        bootstrap_missing = 0;
    otherwise
        error('%s is not a recognized bootstrap strategy', options.bootStrategy);
end
                                                                                                                        
if smc
    %choose the optimal model using smc
    [optTree, results.idxOptTree] = tuning_SMC(results.champions, A, options.n1, options.n2,...
                                                      options.Alpha, results.bootsamples, bootstrap_missing);
else
    %choose the optimal model using a risk function
    [results.idxOptTree, results.fvalues] = tuning_risk(results.prmvalues, results.bootsamples, A, options);
    
    optTree = reults.champions{results.idxOptTree};
end
optP = results.Ps{results.idxOptTree};
    
end
