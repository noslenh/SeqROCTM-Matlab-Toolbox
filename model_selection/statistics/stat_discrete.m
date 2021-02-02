function S = stat_discrete(Br, statistic, threshold)
%STAT_DISCRETE  Test if the branch should be pruned using different
%               statistical criteria
% Inputs 
%
%   Br               : branch to be tested 
%   statistic        : type of statistics used in the pruning criteria. It can
%                        take the values 'context' or 'emp_distribution'
%   threshold        : threshold used in the context algorithm or in the
%                        comparison of the empirical distributions in the 
%
% Outputs
%
%   S   			 : true (prune), H0 is not rejected 
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 05/2019
 
 S = true;
 
 % the leaves from Br that appears less than 2 times won´t be taken into
 % account to do the statistical test
 totals = Br{2,1}';         
 Count = Br{3,1}';          % the transpose is just to coincide with Count in Aline' source code    
 
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
     
     P = bsxfun(@rdivide, Count, totals);
     
     switch statistic
         
         case 'emp_distribution'
             % Compute the distance between the empirical distributions
             for u = 1 : length(P)-1
                 v = u + 1;
                 while (S == true) && (v < length(P) + 1)
                     S = ( abs(P(u)-P(v)) < threshold ) ;
                     v = v + 1;
                 end
                 if S == false
                     break
                 end
             end
             
         case 'context'
             % Compute the statistic given by BIC criteria
             Count_father = sum(Count)'; 
             P_father = Count_father / sum(Count_father);
             
             [~, c] = find(P > 0);   %% inefficient, think a way to avoid the second line
             ind = P > 0;
             
             % based on chi-square asymptotic
             % ss = chi2inv(1-erro, length(Alphabet)-1)/2;
             
             % a particular threshold based on BIC comparison
             % ss = threshold * (length_Alphabet-1) * (d-1) * log(length_X);
             
             ss = threshold;
             test = 2 * sum( Count(ind) .* (log(P(ind)) - log(P_father(c)) ) );
             S =  test < ss ;
     end
 end
            
 end