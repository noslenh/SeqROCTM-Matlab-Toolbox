function d = dissimilarity_emp(ctx, P, est_ctx, est_P, X)
%DISSIMILARITY_EMP Compute a dissimilarity between two probabilistic context tree.
%                 This function compute a dissimilarity between the CTM
%                 (ctx,P) and (est_ctx,est_P) using an empirical estimation
%                 of the invariant measure.
%
% Input
%
% 	ctx           : contexts of the first probabilistic context tree
% 	P             : distributions associated to ctx
% 	est_ctx       : contexts of the second probabilistic context tree
% 	est_P         : distributions associated to est_ctx
% 	X      		  : invariant measure for the context in ctx
%
% Output
%
% 	d             : dissimilarity value
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
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
                fv = get_occurrences(est_ctx{list_v(j)}, X);
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
%GET_OCCURRENCES Get the number of ocurrences of a subsequence in a sequence
%
% Input
%
%   s         : sub-sequence
%   X         : sequence
%
% Output
%
%   f         : number of occurrences
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 04/2019

f = 0;

ls = length(s);
for n = ls : length(X)
    if sum(X(n-ls+1:n) ~= s) == 0
        f = f + 1;
    end
end
end

function p = get_prob(x, contexts, P, mu_ctx)
%GET_PROB Compute the invariant measure of the sub-sequence x  
%
% Input
%
%   x         : sub-sequence to compute the invariant measure
%   contexts  : contexts of the probabilistic context tree
%   P         : distributions associated to contexts
%   mu_ctx    : invariant measure of contexts
%
% Output
%
%   p         : invariant measure for x
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
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
%D_Bhattacharyya Compute the Bhattacharyya distance between two distributions.
%
% Input
%
%   p : first distribution
%   q : second distribution
%
% Output
%
%   d : Bhattacharyya distance
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 02/2019

BC = sum(sqrt(p.*q));
d = - log(BC);

end

function d = D_Hellinger(p, q)
%D_Hellinger Compute the Hellinger distance between two distributions.
%
% Input
%
%   p : first distribution
%   q : second distribution
%
% Output
%
%   d : Hellinger distance
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 04/2019

BC = sum(sqrt(p.*q));
d = sqrt(1 - BC);

end

function d = D_TotalVariation(p, q)
%D_TotalVariation Compute the Total Variation distance between two distributions.
%
% Input
%
%   p : first distribution
%   q : second distribution
%
% Output
%
%   d : Total variation distance
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 04/2019

 d = 0.5 * sum(abs(p - q));

end

function d = D_MatchProbability(p, q)
%D_matchrate Compute the probability that a match happens when symbols are being generated by distributions p and q. 
%
% Input
%
%   p : first distribution
%   q : second distribution
%
% Output
%
%   d : probability of matching
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 08/2021

d = p*q';
end
