function [renwpoint, renewals] = tree_renewalpoint(contexts, P, A, X)
%TREE_RENEWALPOINT Among all the contexts that are renewal points, select the most frequents
%				   and then the smallest. 
%
% Input
%
% contexts 	: set of contexts 
% P 	   	: probability distribution associated to the contexts
% Among	   	: alphabet
% X		   	: sequence of data
%
% Output
%
% renwpoint : renewal point
% renewals  : contexts that are renewal points

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 06/2020

	%
    ncontexts = length(contexts);
    lw = cellfun(@(x) length(x), contexts);
    height = max(lw); 
    
    nA = length(A);
    
	% get the transitions from context to context
    Tw = transition_context2context(contexts, P, A);
    
    %initialize the variable to store the contexts that are renewal points
    renewals = {};
    lrnw = [];
	
	%if for all contexts w, there exist another context v such that v is
    %a suffix of wa, then all the context are renewal points
    ntrans = cellfun(@(x) length(x), Tw);
    multctx = ntrans > 1;
    sw = sum(multctx,2);
    if sum(sw) == 0
        renewals = contexts;
        lrnw = lw;
    else %otherwise, check which contexts are renewal points
        for i = 1 : ncontexts
            if (sw(i) == 0) && (is_renewal(lw(i), i, Tw, height, nA, lw))
                renewals = [renewals, contexts{i}];
                lrnw = [lrnw, lw(i)];
            end
        end
    end
    
    % select the most frequent
    count = countstr(renewals, X);
    freqs = cell2mat(count(1,:));
    [~, idx] = max(freqs);
    
    %choose the smallest
    [~, idxmin]= min(lrnw(idx));
    renwpoint = renewals{idx(idxmin(1))};

end

