function contexts = estimate_functionalSeqROCTM(X, Y, Alphabet, max_height, n_BM, alpha, beta)
%ESTIMATE_FUNCTIONALSEQROCTM estimate a context tree from the sequence of
%                            random objects driven by context tree model (X,Y),
%                            where X is a context tree model and Y is a sequence of functional data
% Inputs
%
% X             : context tree model taking values in A
% Y             : response sequence (sequence of chunk of functions). A matrix of
%                   dimension D x length_X. Each column has a function/vector of
%                   dimension D
% Alphabet      : alphabet 
% max_height    : maximum height of the complete tree
% n_BM          : number of Brownian motion used in the statistical test
% alpha         : significant level of the KS test
% beta          : significance level used in the statistical test with several
%                   Brownian
%
% Outputs
%
% contexts      : estimated context tree

%Author : Noslen Hernandez, Aline Duarte
%Date   : 11/2019

length_X = size(X,2);
[D, Yc] = size(Y);

% validations
if length_X ~= Yc
    error('The length of X and the number of columns of Y must match.')
end

% complete tree
[T, I, nT] = completetree(X, max_height, Alphabet);

if isempty(T)
    contexts = T;
else
         
     % compute the thresholds required in the statitical test used for
     % prunning
     c = sqrt(-1/2 * (log(alpha/2)));
     
     C = binoinv(1-beta, n_BM, alpha);
     
     % project all the functions chunk in the Brownian bridge(s)
     Y_projected = zeros(n_BM, length_X);  % each row contain the projections of all the chunks in a Brownian
     for i = 1 : n_BM
        B = brownianbrigde(D);
        Y_projected(i,:) = dot(Y, B'*ones(1,length_X));
     end

    %
    la = length(Alphabet);
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
    
    for s = max_level-1 : -1 : 1 %iterate the levels buttom-up

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
                
                %
                if isequal(new_node{1,1}, X(1:s-1)) % if X begins with the new_node, the counter needs to be increase by 1
                    new_node{2,1} = new_node{2,1} + 1;
                end
                
                if stat_ks_projective(test(:,b), n_BM, alpha, c, C) == 1 % prune => new_node = leave
                    
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
        % add the phathers of internal nodes as internal nodes for the next
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