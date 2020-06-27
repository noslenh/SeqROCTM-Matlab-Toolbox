function d = dissimilarity_emp(ctx, P, est_ctx, est_P, X)
%DISSIMILARITY_EMP compute a dissimilarity between the two probabilistic context
%              		tree: (ctx,P) and (est_ctx,est_P)
%
% Input
%
% 	ctx           : contexts of the first probabilistic context tree
% 	P             : distributions associated to ctx
% 	est_ctx       : contexts of the second probabilistic context tree
% 	est_P         : distributions associated to est_ctx
% 	A             : Alphabet 
% 	X      		  : invariant measure for the context in ctx
%
% Output
%
% 	d             : dissimilarity value
%
% Usage
%   
%
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 05/2019


%%%%  to-delete %%%%%
est_pi = [];
%%%%%%%%%%%%%%%%%%

d = 0;

emp1 = isempty(ctx);
emp2 = isempty(est_ctx);

if emp1 && emp2     % both trees are empty
    d = D_TotalVariation(P, est_P);
elseif emp1         %ctx is empty and est_ctx not
    for c = 1 : length(est_ctx)
        fv = get_occurrences(est_ctx{c}, X);
        d = d + fv * D_TotalVariation(P, est_P(c,:));
        est_pi = [est_pi, fv];
    end    
elseif emp2         %est_ctx is empty and ctx not
    for c = 1 : length(ctx)
        fv = get_occurrences(ctx{c}, X);
        d = d + fv * D_TotalVariation(P(c,:), est_P);
        est_pi = [est_pi, fv];
    end
else
    match = match_contexts(ctx, est_ctx);
    match = match(:,1);
    
    for c = 1 : length(ctx)
        list_v = match{c};
		lv = length(list_v);
        if lv > 1  % over_estimation
            for j = 1 : lv
                fv = get_occurrences(list_v(j), X);
                d = d + fv * D_TotalVariation(P(c,:), est_P(list_v(j),:));
                est_pi = [est_pi, fv];
            end
        elseif lv == 1 % under_estimation or correct_estimation
            fv = get_occurrences(ctx{c}, X);
            d = d + fv * D_TotalVariation(P(c,:), est_P(list_v,:));
            est_pi = [est_pi, fv];
        end
        % there is a third case lv==0 that only happens if the context does
        % not appear in the sequence, so its fv = 0
    end
end
ss = sum(est_pi);
d = d / ss;
est_pi = est_pi / ss;

end

function f = get_occurrences(s, X)
%    
	f = 0;
	
	ls = length(s);
	for n = ls : length(X)
		if sum(X(n-ls+1:n) ~= s) == 0
			f = f + 1;
		end
	end
end

function p = get_prob(x, contexts, P, mu_ctx)
%GET_PROB compute the invariant measure of the sub-sequence x  
%
% x         : sub-sequence to compute the invariant measure
% contexts  : contexts of the probabilistic context tree
% P         : distributions associated to contexts
% mu_ctx    : invariant measure of contexts
%
% p         : invariant measure for x

%Author : Noslen Hernandez, Aline Duarte
%Date   : 04/2019

    stop = false;
    p = 1;
    l = length(x);
    
    while ~stop
        [w, idx] = contextfunction(x(1:l-1), contexts);
        if ~isempty(w)
            p = p * P(idx, x(l)+1);
            l = l - 1;
        else
            % look for all the contexts that end with x(1:l-1)
            q = 0;
            for c = 1 : length(contexts)
                if (length(contexts{c}) > l-1) && isequal(x(1:l-1), contexts{c}(end-l+2 : end))
                    q = q + P(c, x(l)+1) * mu_ctx(c);
                end
            end
            p = p * q;
            stop = true;
        end
    end
end

function d = D_Bhattacharyya(p, q)
%D_Bhattacharyya compute the Bhattacharyya distance between two distributions
%
% p : first distribution
% q : second distribution
%
% d : Bhattacharyya distance

%Author : Aline Duarte, Noslen Hernandez
%Date   : 02/2019

BC = sum(sqrt(p.*q));
d = - log(BC);

end

function d = D_Hellinger(p, q)
%D_Hellinger compute the Hellinger distance between two distributions
%
% p : first distribution
% q : second distribution
%
% d : Hellinger distance

%Author : Noslen Hernandez, Aline Duarte
%Date   : 04/2019

BC = sum(sqrt(p.*q));
d = sqrt(1 - BC);

end

function d = D_TotalVariation(p, q)
%D_TotalVariation compute the Total Variation distance between two distributions
%
% p : first distribution
% q : second distribution
%
% d : Total variation distance

%Author : Noslen Hernandez, Aline Duarte
%Date   : 04/2019

d = 0.5 * sum(abs(p - q));

end
