%% Example 1: Definition of a context tree model

% alphabet of three symbols
A = [0,1,2];

% context tree with four contexts
tau = {0, 2, [0,1], [1,1]};

% distributions associated to the contexts in tau (4x3 matrix)
% p(0|0)=0, p(1|0)=1, p(2|0)=0 => distribution of context 0
p = [0, 1, 0 ; 1, 0, 0; 0, 0.2, 0.8; 1, 0, 0 ];

% plot the defined context tree
draw_contexttree(tau, A);

%% Example 2: Generating an input sequence using the context tree model defined in Example 1

% length of the sequence
seq_length = 102;

% row vector with the sequence of stimuli (context tree model)
X = generatesampleCTM(tau, p, A, seq_length);

%% Example 3: Load a realization of a SeqROCTM (X,Y) in which the response sequence takes values
%             in the set A, and do model selection using the BIC criteria

% load the input and response sequence (discrete) from file 
load('discrete_data_response.mat');
load('discrete_data_input.mat');

% height of the complete tree
max_height = 6;

% model selection for the finite case using BIC criteria
pen_bic = 1; %penalization constant
[contexts, q] = estimate_discreteSeqROCTM(X, Y, A, max_height,'bic', pen_bic);

% draw the resulting context tree 
figure;
draw_contexttree(contexts, A);

% tuning  the model to choose the BIC constant using the SMC

% values of parameters to estimate the Champion Trees
c_min = 0;
c_max = 30;
[Trees, P, ML, cutoff] = estimate_championTrees2(X, Y, max_height, A, c_min, c_max);

% values of parameters to choose a tree in the Champion Tree using SMC 
n1 = 30; n2 = 90;
alpha = 0.05;
B = 200;

[opt_tree, idtree] = modeltunning_SMC2(Trees, A, 30, 90, alpha, B, 'none', X, [], Trees{1}, P{1}');

% plot a curve of models vs. Likelihood
figure
subplot(1,2,1)
nleaves = cellfun(@(x) size(x,2), Trees);
plot(nleaves, ML, '*--b')
hold on; plot(nleaves(idtree), ML(idtree), 'ro');
text(nleaves(idtree)+0.5, ML(idtree),['\leftarrow C = ' num2str(cutoff(idtree))]);
ylabel('log-likelihood');
xlabel('no. of contexts');

% draw the optimal Tree
subplot(1,2,2)
draw_contexttree(Trees{idtree}, A, [1 0 0], 6);
title('Choosen Model')
