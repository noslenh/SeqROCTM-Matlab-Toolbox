function contexts = estimate_functionalSeqROCTM(X, Y, A, max_height, n_BM, alpha, beta)
%ESTIMATE_FUNCTIONALSEQROCTM Estimate a context tree from a SeqROCTM with functional data as elements of the response sequence.
%           
% Inputs
%
%   X             : context tree model taking values in A
%   Y             : response sequence (sequence of chunk of functions). A matrix of
%                   dimension D x length_X. Each column has a function/vector of
%                   dimension D
%   A             : alphabet 
%   max_height    : maximum height of the complete tree
%   n_BM          : number of Brownian bridges used in the statistical test
%                   or a matrix with the Brownian bridges
%   alpha         : significance level of the KS test
%   beta          : significance level of the Binomial approximation
%
% Outputs
%
%   contexts      : estimated context tree
%
%   References:
%      [1] A. Duarte et al., Mathematic 7, 5 (2019). 
%      [2] N. Hernández et al., arXiv:2009.06371, (2021).
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 01/2021

length_X = size(X,2);
[D, Yc] = size(Y);

% validations
if length_X ~= Yc
    error('The length of X and the number of columns of Y must match.')
end

% complete tree
[T, I] = completetree(X, max_height, A);

if isempty(T)
    contexts = T;
else
     
     if isscalar(n_BM) 
         if n_BM == 0
             B = 1;
             n_BM = size(Y,1);
         else
             % generate the Brownian bridges
             B = zeros(n_BM, D);
             for i = 1 : n_BM
                 B(i,:) = brownianbrigde(D);
             end
         end
     else
         % get the the Brownian bridges 
         B = n_BM;
         n_BM = size(B, 1);
     end
     
     % project all the functions chunk in the Brownian bridge(s)
     % each row contain the projections of all the chunks in a Brownian
     Y_projected = B * Y;
     
     % compute the threshold required in the statistical test used for
     % pruning
     C = binoinv(1-beta, n_BM, alpha);
     
    %
    br_not_test = {};
    
    TEST = cell(max_height+1,1);
    max_level = 0;
    
    % put the complete tree T in an structure TEST organized by levels and
    % by branch (to speed-up)
     for i = 1 : length(T)
        
        lt = length(T{i});
        level = lt + 1;
        if level > max_level, max_level = level; end
        
        % take the projections associated to the context
        prj_w = Y_projected(:,I{i}+lt-1);
        Nw = length(I{i}); 
        
        nr = size(TEST{level},2);
        
        found = false;
        n = 1;
        while ~found && n <= nr
            if isequal(T{i}(2:end), TEST{level}{1,n}{1,1}(2:end))
                found = true;
                TEST{level}{1,n} = [TEST{level}{1,n}, T{i}];    % add the leave
                TEST{level}{2,n} = [TEST{level}{2,n}, Nw];      % add the counts
                TEST{level}{3,n} = [TEST{level}{3,n}, {prj_w}]; % add the associated projections
            else
                n = n + 1;
            end
        end
        if ~found 
            TEST{level}{1, nr + 1} = T(i);
            TEST{level}{2, nr + 1} = Nw;
            TEST{level}{3, nr + 1} = {prj_w};
        end
    end
        
    %
    test = TEST{max_level};
    internal_nodes = {};
    
    for s = max_level-1 : -1 : 1 %iterate the levels bottom-up

        % initialize br_test
        br_test = {};
        for b = 1 : size(TEST{s},2) 
            found = sibling_in_internal_nodes(internal_nodes, TEST{s}{1,b}{1}, 0);
            if found    %not need to test, so add new_node to br_not_test
                br_not_test = [br_not_test, TEST{s}(:,b)];
            else
                br_test = [br_test, TEST{s}(:,b)];
            end
        end
        
        % do statistical test on each branch
        for b = 1 : size(test,2)
             
                % initialize the new_node
                new_node{1,1} = test{1,b}{1}(2:end);
                new_node{2,1} = sum(test{2,b});
                new_node{3,1} = {cell2mat(test{3,b})};
                
                if stat_ks_projective(test(:,b), n_BM, alpha, C) == 1 % prune => new_node = leave
                    
                    % Find if new_node have sibling in the list internal_nodes
                    found = sibling_in_internal_nodes(internal_nodes, new_node{1,1}, 0);
                    
                    if found % not need to test, so add new_node to br_not_test
                        br_not_test = [br_not_test, new_node];
                    else
                        br_test = add_node_to_br_test(br_test, new_node);
                    end
                else % not prune => new_node is an internal_node
                    % add new_node to internal_node list only if a sibling does not
                    % already exist in the list
                    [have_sibling, internal_nodes] = sibling_in_internal_nodes(internal_nodes, new_node{1,1}, 1);
                    if ~have_sibling
                        % if there was no sibling, delete the branch of other sibling
                        % from br_test (if exist) and put in br_not_test.
                        [br, br_test] = delete_branch_from_br_test(br_test, new_node{1,1});
                        br_not_test = [br_not_test, br];
                    end
                    % add the branch already tested to br_not_test
                    br_not_test = [br_not_test, test(:,b)];
                end
        end
        test = br_test;
        % add the fathers of internal nodes as internal nodes for the next
        % iteration
        internal_nodes = cellfun(@(x) x(2:end), internal_nodes, 'UniformOutput', false);
        
    end
    if isempty(br_not_test)
        contexts = {};
    else
        contexts = [br_not_test{1,:}];
    end
end

end

function [have_sibling, internal_nodes] = sibling_in_internal_nodes(internal_nodes, str_node, flag)
% check if a sibling of [str_node] exist in the [internal_nodes] list.
% have_sibling = true if there exist siblings, otherwise have_sibling=false.

% If flag = 1: [str_node] is added to internal_nodes when have_sibling=false

have_sibling = false;
n = 1;

while ~have_sibling && n <= length(internal_nodes)
    if isequal(str_node(2:end), internal_nodes{n}(2:end))
        have_sibling = true;
    else
        n = n + 1;
    end
end

if ~have_sibling && flag, internal_nodes = [internal_nodes; str_node]; end

end

function [br, br_test] = delete_branch_from_br_test(br_test, str_node)
% check if there exist a branch in br_test with siblings of [str_node]. If
% the branch exists, then it is deleted from [br_test] and returned in [br].
% Otherwise, br = [];

found = false;
n = 1;
br = [];
while ~found && n <= size(br_test,2)
    if isequal(br_test{1,n}{1}(2:end), str_node(2:end))
        found = true;
        br = br_test(:,n);
        br_test(:,n) = [];
    else
        n = n + 1;
    end
end
end

function br_test = add_node_to_br_test(br_test, node)
% add [node] to the corresponding branch in [br_test]. If there is no branch
% of sibling, create a new branch with [node].
    
nr = size(br_test,2);
found = false;
n = 1;

while ~found && n <= nr
    if isequal(br_test{1,n}{1,1}(2:end), node{1,1}(2:end))
        found = true;
        br_test{1,n} = [br_test{1,n}, node{1,1}];
        br_test{2,n} = [br_test{2,n}, node{2,1}];
        br_test{3,n} = [br_test{3,n}, node{3,1}];
    else
        n = n + 1;
    end
end

if ~found
    br_test{1, nr + 1} = node(1,1);
    br_test{2, nr + 1} = node{2,1};
    br_test{3, nr + 1} = node{3,1}; 
end

end