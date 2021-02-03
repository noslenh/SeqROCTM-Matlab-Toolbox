% Illustrative example presented in the article XXX

% number of volunteers
n_volunteers = 3;

% probabilistic context tree defining the stimuli
A = [0,1,2];
tau = {[0,0], [1,0], [2,0], [0,1], [1 1], [2,1], 2};
p = [0, 0, 1 ; 0, 0, 1; 0.2, 0.8, 0; 0, 0, 1; 0, 0, 1; 0.2, 0.8, 0; 0.2, 0.8, 0];

% length of the sequences of stimuli
seq_length = 700;

% matrix X of 3x700 containing on each row a sequence of stimuli
Xdata = zeros(3,700);
for v = 1 : n_volunteers
    Xdata(v,:) = generatesampleCTM(tau, p, A, seq_length);
end


% load EEG data for each volunteer
names_volunteer = {'V02', 'V09', 'V19'};

X = [];
Y = cell(1,3);

for v = 1 : n_volunteers
    
    % load stimuli data 
    vname_i = [names_volunteer{v} '_stimuli'];
    x = load(vname_i);
    x = x.(vname_i);
    X = [X; x];

    % load response data
    vname_r = [names_volunteer{v} '_response'];
    y = load(vname_r);
    y = y.(vname_r);
    Y{v} = y;
end

% visualize some symbols of the stimuli sequence and the corresponding EEG
% chunks for volunteer V05
figure;
id_cols = 88:96;
for i = 1 : 9
    % plot the stimuli
    ax = subplot(2, 9, i);
    text(0.5, 0.5, num2str(X(1, id_cols(i))), 'FontSize', 20);
    set( ax, 'visible', 'off')
    
    % plot the EEG chunk
    ax = subplot(2, 9, 9+i);
    plot(Y{1}(:, id_cols(i)));
    set( ax, 'visible', 'off')
    xlim([0 115])
end

% model selection algorithm on the data of each volunteer
nBM = 1000;
alpha = 0.05;
beta = 0.05;

rng(1); tree_v02 = estimate_functionalSeqROCTM(X(1,:), Y{1}, A, 3, nBM, alpha, beta);
rng(1); tree_v09 = estimate_functionalSeqROCTM(X(2,:), Y{2}, A, 3, nBM, alpha, beta);
rng(1); tree_v19 = estimate_functionalSeqROCTM(X(3,:), Y{3}, A, 3, nBM, alpha, beta);

% draw the results
figure
subplot(1,3,1)
draw_contexttree(tree_v02, A, [1 0 0], 3);
subplot(1,3,2)
draw_contexttree(tree_v09, A, [0 1 0], 3);
subplot(1,3,3)
draw_contexttree(tree_v19, A, [0 0 1], 3);

