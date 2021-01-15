function [contexts, Q] = estimate_discreteSeqROCTM(X, Y, Alphabet)
%ESTIMATE_DISCRETESEQROCTM estimate a context tree from the stochastic process driven
%                        by a context tree model (X,Y), where X is a context tree
%                        model and Y is a sequence of discrete data
% Inputs
%
% 	X             : sequence of elements taking values in the alphabet
% 	Y             : response sequence (sequence of elements). A vector of the
%                    same dimension that X
% 	Alphabet      : Alphabet 
%	max_height    : maximum height of the complete tree
% 	statistic     : type of statistics used in the pruning criteria. It can
%                    take the values 'context' or 'emp_distribution'
% 	threshold     : threshold used in the context criterion or in the
%                    emp_distribution criterion
%
% Outputs
%
% 	contexts      : estimated context tree
% 	Q             : estimated family of probability distribution
%
% Usage
%			A = [0,1,2];
%
%			ctxs = {0, [0 1], [1 1], 2}
%			P = [0,   1,   0; ...                  
%     			 0,   0.25,   0.75;
%     			 1,   0,   0;
%     			 1,   0,   0 ];
%
%			X = generatesampleCTM(ctxs, P, A, 100);
%			Y = generatesampleCTM(ctxs, P, A, 100);
%			
%			[c, q] = estimate_discreteSeqROCTM(X, Y, A, 4, 'context', 1);


%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 01/2021

%%%%%%%% name-value pairs arguments
% default values
options = struct('EstimationMethod', 'bic', 'MaxTreeHeight', log(length(X)), ...
                    'ParameterValue', 1, 'CompleteTree', -1, 'TestStructure', -1, 'DegreeOfFreedom', 'fix', ...
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
    df = ~strcmpi(options.DegreeOfFreedom, 'fix');
    if isempty(options.BicPrecomputedStats)
        [contexts, Q] = bic_WCT(X, Alphabet, options.MaxTreeHeight, options.ParameterValue, df, Y);
    else
        [contexts, Q] = bic_WCT(X, Alphabet, options.MaxTreeHeight, options.ParameterValue, df, Y, options.BicPrecomputedStats);
    end
elseif any(strcmpi(options.EstimationMethod, {'context','emp_distribution'}))
    if isequal(options.CompleteTree, -1)
        [contexts, Q] = CTestimator(X, Alphabet, options.MaxTreeHeight, options.EstimationMethod, options.ParameterValue, Y);
    else
        if isequal(options.TestStructure, -1)
            [contexts, Q] = CTestimator(X, Alphabet, options.MaxTreeHeight, options.EstimationMethod, options.ParameterValue, Y, options.CompleteTree);
        else
            [contexts, Q] = CTestimator(X, Alphabet, options.MaxTreeHeight, options.EstimationMethod, options.ParameterValue, Y, options.CompleteTree, options.TestStructure);
        end
    end
else
    error('%s is not a recognized parameter value', options.EstimationMethod);
end

% [contexts, Q] = CTestimator(X, Alphabet, max_height, statistic, threshold, Y); 

end
