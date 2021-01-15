function S = stat_ks_projective(Br, n_BM, alpha, C)
%STAT_KS_PROJECTIVE  Test if the branch should be pruned using a KS test on
%                    the projections of the functional data using n_BM
%                    Brownian motion
% Inputs
%
%   Br    : branch to be tested
%   n_BM  : number of Brownian motion used for the test
%   alpha : significant level of the KS test
%   c     : critical value used in the KS test
%   C     : threshold used to compare the number of rejections when several
%           Brownian motion are used
%
% Output
%
%   S     : true (prune), H0 is not rejected 
%
%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 02/2019

 S = true;      % true => prune 
 
 % delete the leaves from Br that appears less than 2 times
 totals = Br{2,1};         
 projs = Br{3,1};
     
 dlt = totals < 2;
 totals(dlt) = [];
 projs(dlt) = [];
 
 %
 d = length(totals);
 
 if d > 1   % if there exist a branch with at least 2 leaves
     
%      % correct the significant level of the test
%      nt = nchoosek(d, 2);
%      alpha = alpha/nt;
    
     rejections = zeros(n_BM, 1);
     for p = 1 : n_BM
         %statistical test for each Brownian
         reject_H0 = false;
         lv_a = 1;
         lv_b = 2;
         while ~reject_H0 && lv_a < d  %if it is not fulfilled for one pair, then reject H0^w
%          while ~reject_H0 && lv_b < 3 %% invento para hacer solo 1 test     
             % KS test for the pair (lv_a,lv_b)
             
             reject_H0 = kstest2(projs{lv_a}(p,:), projs{lv_b}(p,:), 'Alpha', alpha);
             
%              [~, ~, D] = kstest2(projs{lv_a}(p,:), projs{lv_b}(p,:), 'Alpha', alpha);
%              nrm = sqrt( totals(lv_a)*totals(lv_b) / (totals(lv_a)+totals(lv_b)) );
%              reject_H0 = D * nrm > c;
             
             % update the indices
             if lv_b == d
                 lv_a = lv_a + 1;
                 lv_b = lv_a + 1;
             else
                 lv_b = lv_b + 1;
             end
         end  
         rejections(p) = reject_H0;
     end
     % Check the number of rejections
     if n_BM > 1
         S = sum(rejections) <= C;
     else
         S = ~rejections;
     end
 end
 
end