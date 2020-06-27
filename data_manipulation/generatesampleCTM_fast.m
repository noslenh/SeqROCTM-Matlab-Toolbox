function B = generatesampleCTM_fast(contexts, P, A, lengthSeq, nsamples)
%GENERATESAMPLECTM_FAST generates a samples of context tree model according to
%                       the probabilistic context tree defined by contexts and P
% Inputs
%
%   contexts        : set of contexts
%   P               : probability distributions associated to the contexts. Each
%                     row contains the distribution of the corresponding context
%   A               : Alphabet
%   lengthSeq       : length of the sequences to be generated
%   nsamples        : number of sequences 
%
% Outputs
%
%   B               : a matrix containing on each row a sequence
%
% Usage
%	
%			A = [0,1,2];
%
%			ctxs = {0, [0 1], [1 1], 2}
%			P = [0,   1,   0; ...                  
%     			 0,   0.25,   0.75;
%     			 1,   0,   0;
%     			 1,   0,   0 ];
%
%			B = generatesampleCTM_fast(ctxs, P, A, 100, 300);
%	

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 04/2020

% initialization
B = -1*ones(nsamples, lengthSeq);

if isempty(contexts)    % generate an i.i.d sequence if contexts is empty
    for b = 1 : nsamples
        B(b,:) = samplediscretedist(A, P, lengthSeq);
    end
else
    % get the Markov process representation of the context tree model
    [past, ~, Mc, iT] = contextTree_to_FiniteMarkov(contexts, P, A);
    [npast, max_length] = size(past);
    
    for b = 1 : nsamples
        
        seq = -1*ones(1, lengthSeq);
        
        % initialize the sequence with a past choosen at random
        rp = randperm(npast);
        idx_last_past = rp(1);
        init_past = past(idx_last_past, :); 
        seq(1 : max_length) = init_past;
        next_pos = max_length + 1;

        % add the next symbol according to its probs given the past
        while next_pos < lengthSeq + 1
            [next_symbol, idx] = samplediscretedist(A, Mc(idx_last_past,:), 1);
            seq(next_pos) = next_symbol;
            next_pos = next_pos + 1;
            idx_last_past = iT(idx_last_past, idx);
        end
        
        % save the sequence
        B(b,:) = seq;
    end
end
end