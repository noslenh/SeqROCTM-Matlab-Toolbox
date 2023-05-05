function [optTree, optQ, results] = tune_SeqROCTM(X, Y, A, varargin)
%TUNE_CONTEXTTREEMODEL Tune a SeqROCTM estimation algorithm.
%   [OPTTREE, OPTP] = TUNE_SEQROCTM(X, Y, A) tunes a SeqROCTM using the
%   sequences X,Y with values in the alphabet A. The optimal model model
%   (i.e., context tree and its associated family of distributions) is
%   returned in OPTTREE and OPTQ respectively. 
%
%   [OPTTREE, OPTP, RESULTS] = TUNE_SEQROCTM(...) returns a structure with the following fields:
%       'champions'     -- the Champion Trees
%       'Qs'            -- the family of probability distributions corresponding to each champion tree
%       'fvalues'       -- likelihood or risk function values for the champion trees
%       'prmvalues'     -- the parameter values corresponding to each champion tree 
%       'idxOptTree'    -- index of the optimal model in the champion trees
%       'Xbootsamples'  -- bootstrap samples X generated to tune the model
%       'Ybootsamples'  -- bootstrap samples Y generated to tune the model
%
%   [...] = TUNE_SEQROCTM(X,A,'PARAM1',val1,'PARAM2',val2,...)
%   specifies one or more of the following name/value pairs:
%
%       Parameter                Value
%       'TuningMethod'          'smc' to perform the Smaller Maximizer
%                                Criteria or 'risk' to use a risk function.
%                                Default is 'smc'.
%       'EstimationMethod'       'bic' to estimate the SeqROCTM using
%                                Bayesian Information Criteria and tune the
%                                BIC penalization constant. 
%                               'context_cL' to estimate the SeqROCTM
%                                using the Context Algorithm based on
%                                comparison of likelihoods.
%                                'context_empD' to estimate the SeqROCTM
%                                using the Context Algorithm based on
%                                comparison of distributions. 
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
%       'n1'                     Minimum sample length of bootstrap samples
%                                in SMC. Default value is floor(0.3*length(X)).
%       'n2'                     Maximum sample length of bootstrap samples 
%                                in SMC. Default value is floor(0.9*length(X)).
%       'Alpha'                  Statistical significance of the t-test in
%                                SMC. Default value is 0.01.
%       'BootNSamples'           Number of bootstrap samples. Default is 200.
%       'BootStrategy'           Bootstrap procedure to re-sample the
%                                bivariate sequence (X,Y). 'blocks': a
%                                renewal context is used to create
%                                independent blocks of the sequence (X,Y)
%                                and these blocks are used to create the
%                                bootstrap samples. 'parametric': a model is
%                                used to generate the Y bootstrap samples
%                                from the X bootstrap samples. When this
%                                option is chosen, a bootstrap strategy to
%                                generate the X sequences must be
%                                specified. Default value is 'blocks'.
%       'BootXStrategy'          Bootstrap procedure to re-sample the
%                                sequence X. 'parametric': a model
%                                specified in 'BootXParametricModel' is
%                                used to generate the bootstrap samples.
%                                'blocks': a renewal point is used to
%                                create independent blocks and then sample
%                                from this blocks to generate the bootstrap
%                                samples of the sequence X. 'none': no
%                                bootstrap samples are generate for the
%                                sequence X. Default value is 'none'.
%       'BootXRenewalPoint'      A renewal point to be used when the
%                                'BootXStrategy' is 'block'.
%       'BootXModel'             Context tree model used to generate the
%                                bootstrap samples of the sequence X when
%                                'XBootStrategy' is 'parametric'. A cell
%                                array that contains in BootXModel{1} a
%                                context tree and in BootXModel{2} the
%                                associated transition probabilities. 
%       'BootYModel'             Context tree and distributions used to
%                                generate the bootstrap samples of the
%                                sequence Y. A cell array containing in
%                                BootYModel{1} a context tree and in
%                                BootYModel{2} the distributions. Default
%                                value is [] (in this case the largest
%                                model in the Champion Trees is used). 

%   References:
%      [1] A. Galves et al., Ann. Appl. Stat., 6, 1, 186-209 (2012)
%      [2] N. Hernández et al., arXiv:2009.06371, (2021).  
%      [3] P. Buhlmann et al., Ann. Inst. Statist. Math, 52, 1, 287-315 (2000)
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 02/2021

lX = length(X);

% name-value pairs arguments
% default values
options = struct(   'TuningMethod', 'smc',              ...
                    'EstimationMethod', 'bic',          ...  
                    'MaxTreeHeight', floor(log(lX)),    ...
                    'ParameterLowerBound', 0,           ...
                    'ParameterUpperBound', 100,         ...
                    'Tolerance', 10^-5,                 ...
                    'BicDegreeOfFreedom', 'fix',        ...
                    'n1', floor(0.3*lX),                ...                        
                    'n2', floor(0.9*lX),                ...
                    'Alpha', 0.01,                      ...
                    'BootNSamples', 200,                ...
                    'BootStrategy', 'blocks',           ...
                    'BootXStrategy', 'none',            ...
                    'BootXRenewalPoint', [],            ...
                    'BootXModel', [],                   ...
                    'BootYModel', []                    ...
                    );

% acceptable names
optionNames = fieldnames(options);

for pair = reshape(varargin, 2, [])
    inpName = pair{1};
    
    if any(strcmp(inpName, optionNames))
        
        if strcmpi(inpName, 'estimationmethod')
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

% compute the Champion Trees
[results.champions, results.Ps, results.fvalues, results.prmvalues] = estimate_championTrees2(X, Y, A, ...
                                                                'MaxTreeHeight', options.MaxTreeHeight,             ...
                                                                'ParameterLowerBound', options.ParameterLowerBound, ...
                                                                'ParameterUpperBound', options.ParameterUpperBound, ...
                                                                'EstimationMethod', options.EstimationMethod,       ...
                                                                'Tolerance', options.Tolerance,                     ...
                                                                'BicDegreeOfFreedom', options.BicDegreeOfFreedom    ...
                                                                );

switch options.BootStrategy    
    case 'blocks'    
        
        % define the length of the bootstrap samples
        smc = 1;
        if strcmpi(options.TuningMethod, 'smc')
            boot_samples_length = options.n2;
        elseif strcmpi(options.TuningMethod, 'risk')
            smc = 0;
            boot_samples_length = lX+1;
        else
            error('%s is not a recognized tunning method', options.TuningMethod);
        end
        
        % compute a renewal context of results.champions{1} on X
        [~, Nwa] = countctx(results.champions{1}, X, A);
        p = bsxfun(@rdivide, Nwa, sum(Nwa,2));
        x_renw_point = tree_renewalpoint(results.champions{1}, p, A, X);
        % do the bootstrap based on paired-blocks
        [results.Xbootsamples, results.Ybootsamples] = bootstrap_blocks2(X, Y, x_renw_point, boot_samples_length, options.BootNSamples);
    
    case 'parametric'
        
        % initialize the model used to generate the bootstrapped Y sequences
        if isempty(options.BootYModel)
            tau_y = results.champions{1};
            q_y = results.Ps{1};
        else
            tau_y = options.BootYModel{1};
            q_y = options.BootYModel{2};
        end
        ml = max(cellfun(@(x) length(x), tau_y));
        
        %Define the length of the bootstrap samples
        smc = 1;
        if strcmpi(options.TuningMethod, 'smc')
            initial_Xboot_length = options.n2 + ml;
            boot_samples_length = options.n2;
        elseif strcmpi(options.TuningMethod, 'risk')
            smc = 0;
            initial_Xboot_length = lX+1 + ml;
            boot_samples_length = lX+1;
        else
            error('%s is not a recognized tunning method', options.TuningMethod);
        end
        
        % (1) Generate bootstrap samples for X (three options)
        switch options.BootXStrategy
            case 'blocks'
                if isempty(options.BootXRenewalPoint)
                    error('If you choose "blocks" as bootstrap strategy for the sequence X, you must specify a renewal point.');
                end
                renewal_point = options.BootXRenewalPoint;
                tmp_Xbootsamples = bootstrap_blocks(X, renewal_point, initial_Xboot_length, options.BootNSamples);
            case 'parametric'
                if isempty(options.BootXModel)
                    error('If you choose "parametric" as bootstrap strategy for the sequence X, you must specify a model.');
                end
                tmp_Xbootsamples = generatesampleCTM_fast(options.BootXModel{1}, options.BootXModel{2}, A, initial_Xboot_length, options.BootNSamples);
            case 'none'
                if smc
                    if options.n2 + ml <= lX
                        tmp_Xbootsamples = ones(options.BootNSamples,1) * X(:,1:options.n2+ml);
                    else
                        error('The value for "n2" must be lower than %s', num2str(lX - ml));
                    end
                else
                    error('If TuningMethod is "risk" and BootStrategy is "parametric", you must set BootXStrategy to "parametric". Otherwise, use "risk" with BootStrategy equals "blocks".');
                end
            otherwise
                error('%s is not a recognized bootstrap strategy for the sequence X', options.bootXStrategy);
        end
        
        % Generate bootstrap samples for Y.
        % A parametric bootstrap is used using the bootstrap samples for X
        % and the model tau_y, q_y. This function changes the tmp_Xbootsamples.

        results.Ybootsamples = zeros(options.BootNSamples, boot_samples_length);
        results.Xbootsamples = zeros(options.BootNSamples, boot_samples_length);
        
        i = 1;
        j = 1;
        nsamples = options.BootNSamples;
        while i <= nsamples
            try
                [results.Xbootsamples(j,:), results.Ybootsamples(j,:)] = generatesampleYSeqROCTM(tmp_Xbootsamples(i,:), tau_y, q_y, A);
                j = j + 1;
            catch
                disp(['Eliminating bootstrap sample ' num2str(i) ' due to incompatibilities...']);
                options.BootNSamples = options.BootNSamples - 1;
            end
            i = i + 1;
        end
        results.Ybootsamples(j:end,:) = [];
        results.Xbootsamples(j:end,:) = [];
        
    otherwise
        error('%s is not a recognized bootstrap strategy for the sequence X', options.bootXStrategy);
end

% tuning functions
if smc
    %choose the optimal model using smc
    [optTree, results.idxOptTree] = tuning_SMC2(results.champions, A, options.n1, options.n2, options.Alpha, ...
                                                    results.Xbootsamples, results.Ybootsamples);
else
    %choose the optimal model using a risk function
    [results.idxOptTree, results.fvalues] = tuning_risk2(results.prmvalues, results.Xbootsamples, results.Ybootsamples,...
                                                            A, options);
    
    optTree = results.champions{results.idxOptTree};
end
optQ = results.Ps{results.idxOptTree};
    
end
