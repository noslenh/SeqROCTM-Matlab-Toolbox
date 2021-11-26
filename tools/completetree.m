function [T, I, nT] = completetree(X, max_height, alphabet)
%COMPLETETREE  Compute the complete tree of height max_height compatible
%               with the data X
% Inputs
%
%   X           : sequence of symbols taking values in the alphabet
%   max_height  : height of the complete tree
%   alphabet    : alphabet 
%
% Outputs
%
%   T           : complete tree
%   I           : indexes indicating the position of the contexts of the complete
%                  tree in the sequence X
%   nT          : total number of pairs of siblings in the complete tree (useful
%                  when the prune is based on statistical testing)
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
%			[T, I, nT] = completetree(X, 4, A);
% 

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 05/2020


    T = {};
    I = {};
	nT = 0;
    nson = 0;
    
    % this will take into account all possible past occurring in the sequence
    % X including the past associated to step length(X) (this is the reason
    % why it is written length(X)+1)
    for a = alphabet
        [f, id, nt] = is_leaf(a, alphabet, max_height, 2:length(X)+1, X);
        T = [T, f];
        I = [I, id];
		% update number of pairs of sibling
        nT = nT + nt;
        % counting the sons at the first level 
        if ~isempty(f), nson = nson + 1; end
    end
    if nson > 1
        nT = nT + nchoosek(nson, 2);
    end
end
