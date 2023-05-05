function d = dissimilarity_th(ctx, P, est_ctx, est_P, A, varargin)
%DISSIMILARITY_TH compute a dissimilarity between the two probabilistic context
%                 tree: (ctx,P) and (est_ctx,est_P)
% Input
%
% 	ctx           : contexts of the first probabilistic context tree
% 	P             : distributions associated to ctx
% 	est_ctx       : contexts of the second probabilistic context tree
% 	est_P         : distributions associated to est_ctx
% 	A             : Alphabet 
% 	varargin      : invariant measure for the context in ctx (optional)
%
% Output
%
% 	d             : dissimilarity value
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 02/2021

d = 0;
if ~isempty(varargin)
    mu_ctx = varargin{1};
else
    % compute the invariant measure
    [~, ~, mu_ctx] = EntropyRateCT(ctx, P, A);
end

emp1 = isempty(ctx);
emp2 = isempty(est_ctx);

if emp1 && emp2
    d = D_TotalVariation(P, est_P);
elseif emp1
    for c = 1 : length(est_ctx)
        v = est_ctx{c};
        % compute the probability of v
        d = d + prod(P(v)+1);
    end
elseif emp2
    for c = 1 : length(ctx)
        d = d + mu_ctx(c) * D_TotalVariation(P(c,:), est_P);  
    end
else
    match = match_contexts(ctx, est_ctx);
    match = match(:,1);
    
    for c = 1 : length(ctx)
        v = match{c};
        lv = length(match{c});
        if lv == 0
            % for the context c, there is no over_estimation neither under_estimation
			% because it does not appear in the sequence (this can happen when the sequence is too small)
            d = d + mu_ctx(c); 
        elseif lv > 1   % over_estimation
            for j = 1 : lv
                pv = get_prob(est_ctx{j}, ctx, P, mu_ctx);
                d = d + pv * D_TotalVariation(P(c,:), est_P(v(j),:));
            end
        else            % under_estiamtion
            d = d + mu_ctx(c) * D_TotalVariation(P(c,:), est_P(v,:));
        end
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

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
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

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
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

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
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

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 04/2019

d = 0.5 * sum(abs(p - q));

end
