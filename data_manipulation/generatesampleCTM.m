function [seq, sP] = generatesampleCTM(contexts, P, A, lengthSeq, initStrategy)
%GENERATESAMPLECTM generates a sample of context tree model according to
%                  the probabilistic context tree defined by contexts and P
% Inputs
%
%   contexts        : set of contexts
%   P               : probability distributions associated to the contexts. Each
%                     row contains the distribution of the corresponding context
%   A               : Alphabet
%   lengthSeq       : length of the sequence to be generated
%   initStrategy    : strategy to begin generating the chain. This can take the 
%                     values 'max_length_context' or 'any_string'.
%
% Outputs
%
%   seq       : generated sequence
%   sP        : empirical transition probabilities computed on the generated sample
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
%			X = generatesampleCTM(ctxs, P, A, 100);
%	

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 07/2020

% initialization
seq = -1*ones(1, lengthSeq);
sP = zeros(size(P));
                
if isempty(contexts)    % generate an i.i.d sequence if contexts is empty
    seq = samplediscretedist(A, P, lengthSeq);
    sP = hist(seq, A); %counts in same order than the alphabet
else
    % check if the strategy to begin generating the sequence was specified
    if ~exist('initStrategy', 'var')
        initStrategy = 'max_length_context';
    end
    
    lengths = cellfun(@(x) length(x), contexts);
    max_length = max(lengths); 
    
    switch initStrategy
        case 'max_length_context'
            
             % max_length_context: begin with some of the contexts of maximum length
             % randomly chosen
             idx_ctxs_maxlength = find(lengths == max_length);
             % choose randomly
             msize = numel(idx_ctxs_maxlength);
             rp = randperm(msize);
             idx_ctx = idx_ctxs_maxlength(rp(1));
             ctx_max_length = contexts{idx_ctx};
             
             % initialize the sequence 
             seq(1 : max_length) = ctx_max_length;
             next_pos = max_length + 1;
    
            % add the symbol after the first context
            [next_symbol, idx] = samplediscretedist(A, P(idx_ctx,:), 1);
            seq(next_pos) = next_symbol;
            next_pos = next_pos + 1;
            sP(idx_ctx, idx) = sP(idx_ctx, idx) + 1;
               
        case 'any_string'
            % any_string: begin with any context w randomly chosen.
            % Generate max_length-1 further steps. Truncate the
            % beginning of the sequence deleting m symbols (1<=m<l(w)).
            % This guarantee a sequence that begin with "any" string 
            % (not necessary a context of maximum length, even a context)
            
            % If during the generation of the beginning of the sequence
            % there is no context associated to the current past, we need to
            % begin again
            % In the following we generate a sequence of length lw+max_length-1
            
            success = false;
            while ~success  
               
                % choose randomly the first context
                msize = numel(lengths);
                rp = randperm(msize);
                idx_ctx = rp(1);
                ctx_begin = contexts{idx_ctx};
                lw = length(ctx_begin);

                % initialize the sequence 
                tmp_seq = -1*ones(1, lw + max_length - 1);
                sP = zeros(size(P));
                
                tmp_seq(1 : lw) = ctx_begin;
                next_pos = lw + 1;
                
                % add the symbol after the first context
                [next_symbol, idx] = samplediscretedist(A, P(idx_ctx,:), 1);
                tmp_seq(next_pos) = next_symbol;
                next_pos = next_pos + 1;
                sP(idx_ctx, idx) = sP(idx_ctx, idx) + 1;
                
                trouble = false;
                %while there is a context associated to the past, continue.
                while (next_pos < lw + max_length)&&(~trouble)
                    [~, idx_ctx] = contextfunction(tmp_seq(1:next_pos-1), contexts);
                    if idx_ctx ~= -1
                        [next_symbol, idx] = samplediscretedist(A, P(idx_ctx,:), 1);
                        tmp_seq(next_pos) = next_symbol;
                        next_pos = next_pos + 1;
                        sP(idx_ctx, idx) = sP(idx_ctx, idx) + 1;
                    else
                        trouble = true;
                    end
                end
                % if the sequence was successfully generated we are done
                success = ~trouble;   
            end
            % randomly delete part of the beginning of the sequence
            if lw > 1
                rp = randperm(lw-1);
                ndelete = rp(1);
                tmp_seq(1:ndelete) = [];
                next_pos = next_pos - ndelete;
            end
            seq(1:next_pos-1) = tmp_seq;
    end
    
    % looks for the last context and add the next symbol according to the  
    % distribution associated to such context
    while next_pos < lengthSeq + 1
        [~, idx_ctx] = contextfunction(seq(next_pos - max_length : next_pos-1), contexts);
        [next_symbol, idx] = samplediscretedist(A, P(idx_ctx,:), 1);
        seq(next_pos) = next_symbol;
        next_pos = next_pos + 1;
        sP(idx_ctx, idx) = sP(idx_ctx, idx) + 1;
    end
    
end
% normalize probabilities
sP = bsxfun(@rdivide, sP, sum(sP,2));
end