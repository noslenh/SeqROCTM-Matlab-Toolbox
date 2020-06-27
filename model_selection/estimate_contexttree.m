function [contexts, P] = estimate_contexttree(X, Alphabet, max_height, statistic, bic_pen, varargin)
%ESTIMATE_CONTEXTTREE estimate a context tree from the sequence of dicrete 
%                     symbols X
% Inputs
%
% 	X             : sequence of symbols taking values in Alphabet
% 	Alphabet      : Alphabet 
% 	max_height    : maximum height of the complete tree
% 	statistic     : type of statistics used in the prunning criteria. It can
%                    take the values 'bic' or 'emp_distribution'
% 	bic_pen       : penalization constant used in the BIC criteria or
%                    threshold used in the emp_distribution criteria
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
%			[c, p] = estimate_contexttree(X, A, 4, 'bic', 1);
%
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 05/2019

switch length(varargin)
    case 1
        [contexts, P] = CTestimator(X, Alphabet, max_height, statistic, bic_pen, [], varargin{1});
    case 2
        [contexts, P] = CTestimator(X, Alphabet, max_height, statistic, bic_pen, [], varargin{1}, varargin{2});
    case 0
        [contexts, P] = CTestimator(X, Alphabet, max_height, statistic, bic_pen);
end
end
