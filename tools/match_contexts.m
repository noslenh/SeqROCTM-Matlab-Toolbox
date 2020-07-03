function match = match_contexts(ctx, est_ctx)
%MATCH_CONTEXTS match the contexts in ctx with the corresponding
%               well_estimated, over_estimated or under_estimated contexts
%               in est_ctx
%
% ctx       : set of contexts of the true tree
% est_ctx   : set of contexts of the estimated tree
%
% match     : structure with the matching of contexts
 
%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 04/2019
 
 
nv = length(est_ctx);
nw = length(ctx);
 
match = cell(nw,2);
 
for i = 1 : nv
    v = est_ctx{i};
    lv = length(v);
    for j = 1 : nw
        w = ctx{j};
        lw = length(w);
        if (lw <= lv) && isequal(w, v(end-lw+1 : end))  % overestimation
            match{j,1} = [match{j,1}, i];
            match{j,2} = [match{j,2}, lv - lw];
            break;
        end
        if (lv <= lw) && isequal(v, w(end-length(v)+1 : end)) % underestimation
            match{j,1} = [match{j,1}, i];
            match{j,2} = lv - lw;
        end
    end
end

if isempty(est_ctx)    
    match(:,2) = cellfun(@(x)-length(x), ctx, 'UniformOutput', 0);
end

end