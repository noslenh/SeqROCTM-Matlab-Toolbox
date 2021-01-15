%% Example 1: Definition of a context tree model

% alphabet of three symbols
A = [0,1,2];

% context tree with four contexts
tau = {0, 2, [0,1], [1,1]};

% distributions associated to the contexts in tau (4x3 matrix)
% p(0|0)=0, p(1|0)=1, p(2|0)=0 => distribution of context 0
p = [0, 1, 0 ; 1, 0, 0; 0, 0.2, 0.8; 1, 0, 0 ];

% visualize the context tree used to generate the sequence of stimuli
draw_contexttree(tau, A);

%% Example 2: Generating an input sequence using the context tree model defined in Example 1

% length of the sequence
seq_length = 300;

% row vector with the sequence of stimuli (context tree model)
X = generatesampleCTM(tau, p, A, seq_length);

%% Example 3: Simulate the response sequence of a SeqROCTM (X,Y) and do model selection using the BIC criteria

% Simulate the response data Y for the sequence of inputs generated in
% Example 2, assuming three different strategies for an agent.

% Strategy 1: The agent always that see a 0, plays a 1; always that see a
% 1, plays a 2 and always that see a 2, plays a 0
ctx1 = {0, 1, 2};
q1 = [0 1 0; 0 0 1; 1 0 0];
[X1, Y1] = generatesampleYSeqROCTM(X, ctx1, q1, A);

% Strategy 2: the agent learn the contexts of the input sequence and
% always, at each step, choose the most probable symbol that can came after
% the context identified at such step (this is known as probability
% maximization in neuroscience literature)
ctx2 =  {0, 2, [0,1], [1,1]};
q2 = [0, 1, 0 ; 1, 0, 0; 0, 0, 1; 1, 0, 0 ]; 
[X2, Y2] = generatesampleYSeqROCTM(X, ctx2, q2, A);

% strategy 3: The agent does not learn anything and play randomly choosing
% at each step any symbol of the alphabet independently and in a uniform
% way.
ctx3 =  {};
q3 = [1/3 ; 1/3; 1/3 ]; 
[X3, Y3] = generatesampleYSeqROCTM(X, ctx3, q3, A);

% From the data (X,Y) estimate the parameters (ctx,q) that describe the
% strategy of the agent using the model selection algorithms

% parameters value to estimate the Champion Trees
c_min = 0;
c_max = 50;     %if the value is not high enough the function returns a warning
max_height = 6;

% parameters value to select a tree from the Champion Tree set using SMC 
n1 = ceil(0.3*seq_length) ; n2 = ceil(0.9*seq_length);
alpha = 0.05;
B = 200;

% Estimate Champion Trees on each case
[Trees1, P1, ML1, cutoff1] = estimate_championTrees2(X1, Y1, max_height, A, c_min, c_max);
[Trees2, P2, ML2, cutoff2] = estimate_championTrees2(X2, Y2, max_height, A, c_min, c_max);
[Trees3, P3, ML3, cutoff3] = estimate_championTrees2(X3, Y3, max_height, A, c_min, c_max);

% Apply SMC to choose the optimal model on each case
[opt_tree1, idtree1] = modeltunning_SMC2(Trees1, A, n1, n2, alpha, B, 'none', X1, [], Trees1{1}, P1{1}');

[opt_tree2, idtree2] = modeltunning_SMC2(Trees2, A, n1, n2, alpha, B, 'none', X2, [], Trees2{1}, P2{1}');
% %if the context tree model generating the input sequence is known, it can
% %used other bootstrap strategies for the input sequence (this can be useful
% %when the sample size is small), e.g.,
% renewalpoint = tree_renewalpoint(contexts, p, A, X);
% [opt_tree2, idtree2] = modeltunning_SMC2(Trees2, A, n1, n2, alpha, B, 'blocks', X2, renewalpoint, Trees2{1}, P2{1}');

[opt_tree3, idtree3] = modeltunning_SMC2(Trees3, A, n1, n2, alpha, B, 'none', X3, [], Trees3{1}, P3{1}');

% show the results of the estimation procedures
figure
for i = 1 : 3
    subplot(2,3,i)
    eval(['Trees = Trees' num2str(i) ';']);
    eval(['ML = ML' num2str(i) ';']);
    eval(['idtree = idtree' num2str(i) ';']);
    eval(['cutoff = cutoff' num2str(i) ';']);   
    nleaves = cellfun(@(x) size(x,2), Trees);
    plot(nleaves, ML, '*--b')
    hold on; plot(nleaves(idtree), ML(idtree), 'ro');
    text(nleaves(idtree)+0.5, ML(idtree), ['\leftarrow C = ' num2str(cutoff(idtree))], 'FontSize', 8);
    ylabel('log-likelihood');
    xlabel('no. of contexts');
    % draw the choosen context trees
    subplot(2,3,3+i)
    draw_contexttree(Trees{idtree}, A, [1 0 0], 3);
end

% Calling the model selection procedure without tuning the bic constant
% (as we are using the optimal value for the bic constant obtained during
% the tuning, the models will be the same than the obtained before)
[tau1, q1] = estimate_discreteSeqROCTM(X1, Y1, A, max_height, 'bic', cutoff1(idtree1));
[tau2, q2] = estimate_discreteSeqROCTM(X2, Y2, A, max_height, 'bic', cutoff2(idtree2));
[tau3, q4] = estimate_discreteSeqROCTM(X3, Y3, A, max_height, 'bic', cutoff3(idtree3));
