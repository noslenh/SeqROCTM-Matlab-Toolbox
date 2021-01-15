function [contexts, P] = estimate_contexttree(X, Alphabet, varargin)
%ESTIMATE_CONTEXTTREE estimate a context tree from the sequence of discrete 
%                     symbols X
% Inputs
%
% 	X             : sequence of symbols taking values in Alphabet
% 	Alphabet      : Alphabet 
% 	max_height    : maximum height of the complete tree
% 	statistic     : type of statistics used in the pruning criteria. It can
%                    take the values 'context' or 'emp_distribution'
% 	threshold     : threshold used in context criterion or in the
%                    emp_distribution
%   varargin{1}   : complete tree (contexts and indexes)
%   varargin{2}   : TEST structure

%
% Outputs
%
% 	contexts      : estimated context tree
% 	P             : estimated family of probability distribution
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
%			
%			[c, p] = estimate_contexttree(X, A, 'MaxTreeHeight', 4, 'EstimationMethod', 'bic', 'ParameterValue', 2);
%
%

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
    df = ~strcmpi(options.DegreeOfFreedom,'fix');
    if isempty(options.BicPrecomputedStats)
        [contexts, P] = bic_WCT(X, Alphabet, options.MaxTreeHeight, options.ParameterValue, df);
    else
        [contexts, P] = bic_WCT(X, Alphabet, options.MaxTreeHeight, options.ParameterValue, df, [], options.BicPrecomputedStats);
    end
elseif any(strcmpi(options.EstimationMethod, {'context','emp_distribution'}))
    if isequal(options.CompleteTree, -1)
        [contexts, P] = CTestimator(X, Alphabet, options.MaxTreeHeight, options.EstimationMethod, options.ParameterValue);
    else
        if isequal(options.TestStructure, -1)
            [contexts, P] = CTestimator(X, Alphabet, options.MaxTreeHeight, options.EstimationMethod, options.ParameterValue, [], options.CompleteTree);
        else
            [contexts, P] = CTestimator(X, Alphabet, options.MaxTreeHeight, options.EstimationMethod, options.ParameterValue, [], options.CompleteTree, options.TestStructure);
        end
    end
else
    error('%s is not a recognized parameter value', options.EstimationMethod);
end


% switch length(varargin)
%     case 1
%         [contexts, P] = CTestimator(X, Alphabet, max_height, statistic, threshold, [], varargin{1});
%     case 2
%         [contexts, P] = CTestimator(X, Alphabet, max_height, statistic, threshold, [], varargin{1}, varargin{2});
%     case 0
%         [contexts, P] = CTestimator(X, Alphabet, max_height, statistic, threshold);
% end
end
