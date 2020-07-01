function [contexts, P] = CTestimator(X, Alphabet, max_height, statistic, threshold, varargin)
%CTESTIMATOR estimate a context tree from the sequence X or from the 
%            pair of sequences (X,Y) -- sequence of random objects driven by a context tree model
% Inputs
%
% X             : sequence of symbols taking values in Alphabet
% Y             : response sequence (sequence of symbols). A vector of the
%                   same dimension that X (it is optionally given in varargin)
% Alphabet      : Alphabet 
% max_height    : maximum height of the complete tree
% statistic     : type of statistics used in the prunning criteria. It can
%                   take the values 'bic' or 'emp_distribution'
% threshold     : penalization constant used in the BIC criteria or
%                   threshold used in the comparison of the empirical distributions
% varargin{1}   : Y sequence
% varargin{2}   : complete tree (contexts and indexes)
% varargin{3}   : TEST structure
%
% Output
%
% contexts      : estimated context tree
% P             : estimated family of probability distribution
%
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
%			[c, p] = CTestimator(X, A, 4, 'bic', 1, Y);
%
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 05/2019

compute_complete_tree = true;
compute_TEST = true;

% check the values of varargin
switch length(varargin)
    case 1
        Y = varargin{1};
    case 2
        Y = varargin{1};
        T = varargin{2}{1};
        I = varargin{2}{2};
        compute_complete_tree = false;
    case 3
        Y = varargin{1};
        T = varargin{2}{1};
        I = varargin{2}{2};
        TEST = varargin{3};
        max_level = length(TEST);
        compute_complete_tree = false;
        compute_TEST = false;
    case 0
        Y = X;
end

length_X = length(X);

% compute the complete tree
if compute_complete_tree
    [T, I] = completetree(X, max_height, Alphabet);
end

if isempty(Y), Y = X; end

if isempty(T)   % if the complete tree is the empty tree
    contexts = T;
    counts = histc(Y, Alphabet);
    P = counts / sum(counts); 	% return the frequency of each symbol
    P = P';
else
    la = length(Alphabet);
    br_not_test =  {};
    
    % organize the information in a structure called TEST to facilitate the
    % prunning process
    if compute_TEST
        TEST = cell(max_height+1,1);
        max_level = 0;
    
        % for each leaf in the complete tree
        for i = 1 : length(T)
        
            level = length(T{i}) + 1;
            if level > max_level, max_level = level; end
            
            % get the number of times the leaf has ocurred and the
            % transitions
            [Nw, Nwa] = get_counts(T{i}, I{i}, Y, la);
            
            % number of elements in that level
            nr = size(TEST{level},2);
        
            found = false;
            n = 1;
            % looks for some sibling in the level. If found it, put it
            % together
            while ~found && n <= nr
                if isequal(T{i}(2:end), TEST{level}{1,n}{1,1}(2:end))
                    found = true;
                    TEST{level}{1,n} = [TEST{level}{1,n}, T{i}];  % add the leaf
                    TEST{level}{2,n} = [TEST{level}{2,n}, Nw];    % add the counts
                    TEST{level}{3,n} = [TEST{level}{3,n}, Nwa];   % add the transitions  
                else
                    n = n + 1;
                end
            end
            % if no sibling were found, create a new item 
            if ~found 
                TEST{level}{1, nr + 1} = T(i);
                TEST{level}{2, nr + 1} = Nw;
                TEST{level}{3, nr + 1} = Nwa;
            end
        end
    end
        
    % Here begin the prunning procedure
    test = TEST{max_level};
    internal_nodes = {};
    
    for s = max_level-1 : -1 : 1 
        
        % initialize
        flag_leaf_without_occurrences = false;
        
        % initialize the branch to be tested, br_test
        br_test = {};
        
        % for each branch in TEST (i.e., set of siblings), if one of the
        % sibling is an internal not, we do not need to test that branch
        for b = 1 : size(TEST{s},2) 
            found = sibling_in_internal_nodes(internal_nodes, TEST{s}{1,b}{1}, 0);
            if found % not need to test, so add it to br_not_test
                br_not_test = [br_not_test, TEST{s}(:,b)];
            else
                br_test = [br_test, TEST{s}(:,b)];
            end
        end
        
        % test each of the branch to be tested
        for b = 1 : size(test,2)
             
                % do the new_node (father)
                new_node{1,1} = test{1,b}{1}(2:end);
                new_node{2,1} = sum(test{2,b});
                new_node{3,1} = sum(test{3,b},2);
                
                % TRICK: if X begins with new_node, the counters needs to be increased by 1
                % because that instance of new_node it is not taken into account by the sons 
                if isequal(new_node{1,1}, X(1:s-1))
                    % increasse the frequency
                    new_node{2,1} = new_node{2,1} + 1;
                    % index of the next symbol in the alphabet
                    idx = Y(s) + 1; 
                    % increases the transition
                    new_node{3,1}(idx) = new_node{3,1}(idx) + 1;
                end
                
                % call the statistics
                if stat_discrete(test(:,b), statistic, length_X, la, threshold) == 1 % prune => new_node = leave
                    
                    %%%% FIRST, verify if new_node has positive frequency up to lenght(X)-1, i.e., if Nw>0     %%%% 
                    %%%% This is a very low probability event, but it can happen when the sample size is small %%%%
                    %%%% and the unique ocurrence of new_node was at the end of the sequence, so Nw = 0        %%%%
                    if new_node{2,1} == 0
                        flag_leaf_without_occurrences = true;
                        father_lwo = new_node{1,1}(2:end);
                        
                        % we need to prune everything that has branched from the father of new_node
                        % To do this:
                        % 1) delete siblings of new_node from internal_nodes
                        lf = length(father_lwo);
                        to_delete = [];
                        for ii = 1 : length(internal_nodes)
                            if isequal(father_lwo, internal_nodes{ii}(end-lf+1:end))
                                to_delete = [to_delete; ii];
                            end
                        end
                        internal_nodes(to_delete) = [];
                        
                        % 2) delete siblings of new_node from br_test
                        [~, br_test] = delete_branch_from_br_test(br_test, new_node{1,1});
                        
                        % 3) delete from br_not_test every context that has the father of new_node as a suffix
                        to_delete = [];
                        for ii = 1 : size(br_not_test,2)
                            if iscell(br_not_test{1,ii})
                                suff = br_not_test{1,ii}{1}(end-lf+1:end);
                            else
                                suff = br_not_test{1,ii}(end-lf+1:end);
                            end
                            if isequal(father_lwo, suff)
                                to_delete = [to_delete; ii];
                            end
                        end
                        br_not_test(:,to_delete) = [];
                       
                        % 4) add to br_test a node with the statistics of the father of new_node. In the next
                        % iteration this node is pruned by the statistical test (compute its frequency and  
                        % transition from X)
%                         node{1,1} = father_lwo;
%                         node{2,1} = 0;
%                         node{3,1} = zeros(la,1);
                        for ii = 1 : length_X - lf 
                            if isequal(father_lwo, X(ii:ii+lf-1)) % one step before the end of the sequence
                                new_node{2,1} = new_node{2,1} + 1;
                                new_node{3,1}(Y(ii+lf)+1) = new_node{3,1}(Y(ii+lf)+1) + 1;
                            end
                        end
                        br_test = add_node_to_br_test(br_test, new_node);
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % Find if new_node have sibling in internal_nodes
                    found = sibling_in_internal_nodes(internal_nodes, new_node{1,1}, 0);
                    
                    if found % not need to test, so add new_node to br_not_test
                        br_not_test = [br_not_test, new_node];
                        % When flag_leaf_without_occurrences is true, found is always false, except if new_node is
                        % from a different branch. In such a case, the procedure is as ussual.
                    else
                        % 
                        if (~flag_leaf_without_occurrences)||(~isequal(father_lwo, new_node{1,1}(2:end)))
                            br_test = add_node_to_br_test(br_test, new_node);
                        end
                    end
                else % not prune => new_node is an internal_node
                    % add new_node in internal_node only if a sibling does not exist already in internal_node
                    [have_sibling, internal_nodes] = sibling_in_internal_nodes(internal_nodes, new_node{1,1}, 1);
                    if ~have_sibling
                        if (~flag_leaf_without_occurrences)||(~isequal(father_lwo, new_node{1,1}(2:end)))
                            % if there was no sibling in internal nodes, delete the branch of other sibling
                            % from br_test (if exist) and add it to br_not_test (contexts)
                            [br, br_test] = delete_branch_from_br_test(br_test, new_node{1,1});
                            br_not_test = [br_not_test, br];
                            % add the branch already tested to br_not_test
                            br_not_test = [br_not_test, test(:,b)];
                        else
                            % delete new_node from internal_nodes (it is the last one, because it was the last addition)
                            internal_nodes{end} = [];
                        end
                    else
                        % Here flag_leaf_without_occurrences is always false, except if new_node is
                        % from a different branch. In such a case, the procedure is as ussual.
                        
                        % add the branch already tested to br_not_test
                        br_not_test = [br_not_test, test(:,b)];
                    end  
                end
        end
        test = br_test;
        % add father of internal nodes as internal nodes for the next
        % iteration
        internal_nodes = cellfun(@(x) x(2:end), internal_nodes, 'UniformOutput', false);
        
    end
    if isempty(br_not_test)
        contexts = {};
        counts = histc(Y, Alphabet);
        P = counts / sum(counts); % return the frequency of each symbol
        P = P';
    else
        contexts = [br_not_test{1,:}];
        P = [br_not_test{3,:}];
        P = bsxfun(@rdivide, P, [br_not_test{2,:}]);
    end
end

end

function [have_sibling, internal_nodes] = sibling_in_internal_nodes(internal_nodes, str_node, flag)
% check if a sibling of [str_node] exist in the [internal_nodes] list.
% have_sibling = true if it exists, otherwise have_sibling = false.

% If flag = 1: node is added to internal_nodes when have_sibling = false

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
% check if there exist a branch in br_test with siblings of [str_node]. If it
% exists, then such brach is deleted from [br_test] and return in [br].
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
% add a node to the correspondind branch in br_test. If there is no branch
% of sibling, create a new branch with node
    
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

function [Nw, Nwa] = get_counts(w, ind, X, length_alphabet)
%GET_COUNTS Gives how often w appears in X and how often each symbol of the
%           alphabet appears after w.
%
% Inputs
%   w               : sub-sequence
%   ind             : position where w happens in the sequence
%   X               : sequence
%   length_alphabet : length of the alphabet
%
% Outputs
%   Nw              : Number of ocurrences of w in X
%   Nwa             : Number of occurences of each symbol in the alphabet
%                     after w 

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 05/2019

    Nwa = zeros(length_alphabet,1);
    lw = length(w);
    lx = length(X);
    
    for i = 1 : length(ind)
        pos = ind(i) + lw;
        if pos <= lx            % this is because ind+l(w) gives the position after w
            loc = X(pos) + 1;   % faster way: interpreting the symbol as index
            Nwa(loc) = Nwa(loc) + 1;
        end  
    end
    
    Nw = sum(Nwa);

end
