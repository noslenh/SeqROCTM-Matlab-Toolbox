function [contexts, Q] = estimate_discreteSeqROCTM( X, Y, Alphabet, max_height, statistic, bic_pen)
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
% 	statistic     : type of statistics used in the prunning criteria. It can
%                    take the values 'bic' or 'emp_distribution'
% 	bic_pen       : penalization constant used in the BIC criteria or
%                    threshold in the emp_distribution criteria
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
%			[c, q] = estimate_discreteSeqROCTM(X, Y, A, 4, 'bic', 1);


%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 05/2019

[contexts, Q] = CTestimator(X, Alphabet, max_height, statistic, bic_pen, Y); 

end
