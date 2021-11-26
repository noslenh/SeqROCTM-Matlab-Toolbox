function [X, Y, qemp] = generatesampleXYSeqROCTM(ctx_X, p, ctx_Y, q, A, n)
%GENERATESAMPLEYSEQROCTM Generates a SeqROCTM using a CTM to generate the
%                        input sequence and a context tree with the
%                        distributions to generate the response sequence
%
% Inputs
%
% ctx_X    : set of contexts to generate the input sequence X
% p        : distributions associated to contexts in ctx_X
% ctx_Y    : set of contexts to generate the response sequence Y
% q        : distributions associated to contexts in ctx_Y
% A        : alphabet
% n        : length of the sequences to be generated
%
% Outputs
%
% X        : sequence of inputs
% Y        : sequence of responses
% qemp     : empirical distributions computed on the simulated samples
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 06/2020
    
    % generate the input sequence
    max_lengthX = max(cellfun(@(x) length(x), ctx_X));
    X = generatesampleCTM(ctx_X, p, A, n + max_lengthX, 'max_length_context');
    
    % generate the response sequence
    [X, Y, qemp] = generatesampleYSeqROCTM(X, ctx_Y, q, A);
    
end