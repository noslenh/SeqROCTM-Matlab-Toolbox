function [Trees, P, ML, cutoff] = estimate_championTrees(X, max_height, A, varargin) %criterion
%ESTIMATE_CHAMPIONTREES estimate a context tree using the 'bic' criterion
%                       for different values of the penalization constant.
%                       The set of trees retrieved for all value of the
%                       penalization constant is called champion trees (see
%                       xxx for more details)
%
% Input
%   
%   X           : sequence of data
%   max_height  : height of the complete tree
%   A           : alphabet
%   varargin    : contain the l_min and u values
%                 l_min: minimum value for the penalization constant (usually 0)
%                 u    : maximum value for the penalization constant
%                        (usually a big value such that the resulting tree is the empty tree)
%                 tol  : tolerance used when estimating the set of champion trees  
%
% Output
%   
%   Trees       : set of champion trees
%   P           : set of distributions associated to the trees
%   ML          : maximum likelihood value for each of the champion tree
%   cuttoff     : values of the bic penalization
%
% Usage
%
%

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 05/2020


numvarargs = length(varargin);
optargs = {0, 100, 10^-5};
optargs(1:numvarargs) = varargin;            
[l_min, u, tol] = optargs{:};

% compute the complete tree and the TEST structure only once (for speed-up)
[T, I] = completetree(X, max_height, A);
TEST = getTESTstructure(T, I, max_height, length(A), X);

% estimate the trees for the minimum and maximal value of the penalization
% constant
[tau_l, p_l] = estimate_contexttree(X, A, max_height, 'bic', l_min, {T, I}, TEST);  
[tau_upper, p_upper] = estimate_contexttree(X, A, max_height, 'bic', u, {T, I}, TEST);

% initialize
upper_bound = u;
Trees = {}; P = {}; ML = []; cutoff = [];

if ~isempty(tau_upper)
    disp('Warning: The empty tree is not obtain for the maximun value of the penalization constant given.')
end

    tau_u = tau_upper;
    p_u = p_upper;

    i = 1; Trees{i} = tau_l; P{i} = p_l; ML(i) = treeloglikelihood(tau_l, X, A); cutoff(i) = l_min;

    % estimate the different trees in the interval specified for the
    % penalization constant
    while ~isequalCT(tau_l, tau_upper)
        while abs(u - l_min) > tol
            while ~isequalCT(tau_u, tau_l)&&(abs(u - l_min) > 10^-5) % the second condition its necessary because for some cases,
                a = u;                                               % the complete tree is obtain when l_min=0 and for any value     
                tau_a = tau_u; p_a = p_u;                            % greater than zero, a tree different from the complete tree is obtained   
                u = (l_min + u)/2;                                      
                [tau_u, p_u] = estimate_contexttree(X, A, max_height, 'bic', u, {T, I}, TEST);   
            end
            l_min = u; tau_l = tau_u;
            u = a; tau_u = tau_a; p_u = p_a;
%             [tau_u, p_u] = estimate_contexttree(X, A, max_height, 'bic', u, {T, I}, TEST);
        end
        i = i + 1;
        Trees{i} = tau_u; P{i} = p_u; cutoff(i) = u;
        ML(i) = treeloglikelihood(tau_u, X, A);
        l_min = u; tau_l = tau_u; 
        u = upper_bound;
        tau_u = tau_upper;  
        p_u = p_upper;
    end
end