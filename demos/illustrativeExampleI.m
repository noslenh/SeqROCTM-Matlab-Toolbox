%% Definition of a context tree model

% alphabet of three symbols
A = [0,1,2];

% context tree containing four contexts 0, 2, 01, 11
tau = {0, 2, [0,1], [1,1]};

% distributions associated to each contexts in tau (4x3 matrix)
% p(0|0)=0, p(1|0)=1, p(2|0)=0 => distribution of context 0
p = [0, 1, 0 ; 1, 0, 0; 0, 0.2, 0.8; 1, 0, 0 ];

% visualize the context tree used to generate the sequence of stimuli
draw_contexttree(tau, A);

%% Generation of an input sequence using the context tree model defined above

% length of the sequence
seq_length = 300;

% row vector with the sequence of stimuli (context tree model)
X = generatesampleCTM(tau, p, A, seq_length);

%% Simulation of the response sequence of a SeqROCTM (X,Y) 

% Simulate the response data Y for the sequence of inputs generated in
% Example 2, assuming three different strategies for an agent.

% Strategy 1: The agent always that see a 0, plays a 1;
%             the agent always that see a 1, plays a 2; and
%             the agent always that see a 2, plays a 0.

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

% Strategy 3: The agent does not learn anything and plays randomly choosing
% at each step any symbol from the alphabet independently and in a uniform
% way.
ctx3 =  {};
q3 = [1/3 ; 1/3; 1/3 ]; 
[X3, Y3] = generatesampleYSeqROCTM(X, ctx3, q3, A);

%%  Model selection

% From the data (X,Y) estimate the parameters (ctx,q) that describe the
% strategy of the agent using the model selection algorithms

% some parameters value
c_min = 0;
c_max = 1000;     %if the value is not high enough, the function returns a warning message
max_height = 6;
 
alpha = 0.05;

% tune the SeqROCTM model for each strategy
[~,~, r1] = tune_SeqROCTM(X1, Y1, A, 'TuningMethod', 'smc',             ...
                                     'EstimationMethod', 'context_empD',...
                                     'MaxTreeHeight', max_height,       ...
                                     'ParameterLowerBound', c_min,      ...
                                     'ParameterUpperBound', c_max,      ...
                                     'Alpha', alpha);
                                 
[~,~, r2] = tune_SeqROCTM(X2, Y2, A, 'TuningMethod', 'smc',             ...
                                     'MaxTreeHeight', max_height,       ...
                                     'EstimationMethod', 'context_cL',  ...
                                     'ParameterLowerBound', c_min,      ...
                                     'ParameterUpperBound', c_max,      ...
                                     'Alpha', alpha);
                                 
[~,~, r3] = tune_SeqROCTM(X3, Y3, A, 'TuningMethod', 'smc',             ...
                                     'MaxTreeHeight', max_height,       ...
                                     'EstimationMethod', 'bic',         ...
                                     'ParameterLowerBound', c_min,      ...
                                     'ParameterUpperBound', c_max,      ...
                                     'Alpha', alpha,                    ...
                                     'BootNSamples', 200,               ...
                                     'BootStrategy', 'blocks');

% show the results of the estimation procedures
figure
for i = 1 : 3
    subplot(2,3,i)
    % get the structure of the corresponding model
    eval(['r = r' num2str(i) ';']); 
    % get the values from the structure r
    nleaves = cellfun(@(x) size(x,2), r.champions);
    ML = r.fvalues;
    idtree = r.idxOptTree;
    cutoff = r.prmvalues;
    % draw the curve
    plot(nleaves, ML, '*--b')
    hold on; plot(nleaves(idtree), ML(idtree), 'ro');
    text(nleaves(idtree)+0.5, ML(idtree), ['\leftarrow C = ' num2str(cutoff(idtree))], 'FontSize', 8);
    ylabel('$$\log(L_{(\tau, \hat{q})}(Y_1^n|X_1^n))$$', 'interpreter', 'latex');
    xlabel('$$|\tau|$$', 'interpreter', 'latex');
    % draw the chosen context trees
    subplot(2,3,3+i)
    draw_contexttree(r.champions{idtree}, A, [1 0 0], 3);
end

% Calling the model selection procedure without tuning (using the default
% value of the hyper-parameter)
[tau1, q1] = estimate_discreteSeqROCTM(X1, Y1, A, 'MaxTreeHeight', max_height, 'EstimationMethod', 'context_empD');
[tau2, q2] = estimate_discreteSeqROCTM(X2, Y2, A, 'MaxTreeHeight', max_height, 'EstimationMethod', 'context_cL');
[tau3, q3] = estimate_discreteSeqROCTM(X3, Y3, A, 'MaxTreeHeight', max_height, 'EstimationMethod', 'bic', 'ParameterValue', 0.08);

% show the results in the console
print_tree(tau1);
print_tree(tau2);
print_tree(tau3);
