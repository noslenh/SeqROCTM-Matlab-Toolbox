function S = stat_discrete(Br, statistic, threshold)
%STAT_DISCRETE  Test if the branch should be pruned using different statistical criteria
% Inputs 
%
%   Br               : branch to be tested 
%   statistic        : type of statistics used in the pruning criteria. It can
%                        take the values 'context_cL' or 'context_empD'
%   threshold        : threshold used in the context algorithm or in the
%                        comparison of the empirical distributions in the 
%
% Outputs
%
%   S   			 : true (prune), H0 is not rejected 
%
%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 12/2022
 
 S = true;
 
 % the leaves from Br that appears less than 2 times will not be taken into
 % account to do the statistical test
 totals = Br{2,1}';         
 Count = Br{3,1}';              
 
 lt = length(totals);
 dlt = false(lt,1);
 leaf_data = true;
 for i = 1 : lt
     if totals(i) < 2   % THIS FILTER THE LEAVES TAKEN INTO ACCOUNT
         dlt(i) = true;
         if totals(i) == 0
             leaf_data = false; % if at least one leaf have no data, prune the branch
         end
     end
 end

 totals(dlt) = [];
 Count(dlt,:) = [];
 
 % number of leaves to do the statistical test
 d = length(totals);
 
 if (d > 1)&&(leaf_data)
     
     % distributions associated to the leaves
     P = bsxfun(@rdivide, Count, totals);
     
     % distribution associate to the father
     Count_father = sum(Count)'; 
     P_father = Count_father / sum(Count_father);
     
     switch statistic
         
         case 'context_empD'
             
             % Compute the distance between the empirical distributions
             S =  max ( max( abs(P - ones(d,1) * P_father') ) ) < threshold;
             
         case 'context_cL'
             
             % Compute the statistic using likelihoods
             [~, c] = find(P > 0);   %% inefficient, think a way to avoid the second line
             ind = P > 0;
             
             ss = threshold;
             test = 2 * sum( Count(ind) .* (log(P(ind)) - log(P_father(c)) ) );
             S =  test < ss ;
             
             % based on chi-square asymptotic
             % ss = chi2inv(1-erro, length(Alphabet)-1)/2;
             
             % a particular threshold based on BIC comparison
             % ss = threshold * (length_Alphabet-1) * (d-1) * log(length_X);
     end
 end
            
 end