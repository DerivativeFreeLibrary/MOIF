% DMS Version 0.2.
%
% Copyright (C) 2011 A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz,
% and L. N. Vicente.
%
% http://www.mat.uc.pt/dms

% MOIF v.0.1
% Copyright (C) 2017 G.Cocchi, G.Liuzzi, A.Papini, M.Sciandrone:


% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
error_params = false;

% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.% Output Options.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
output = 1; % 0-1 variable: 0 if only a final report is displayed at the
% screen; 1 if at each iteration output is displayed at the
% screen and recorded in a text file stored at the current directory.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stopping Criteria.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
purityPercentage = 1; % Purity percentage coefficient (see paper). It has to be a number between 0 and 1.

if ~check_params(purityPercentage, 'double', true)
    disp('Purity percentage is a double non-negative param')
    
end
% If 1 a line search step is always tried, otherwise never.

box_projections = false; % Projections on box constraints before evaluating the stencil.

if ~check_params(box_projections, 'logical', false)
    disp('Box_projection percentage is a logical param')
    error_params = true;
end


tol_stop  = 10^-3; % Minimum stepsize for stencil evaluation.

if ~check_params(tol_stop, 'double', true)
    disp('tol_stop is a double non-negative param')
    error_params = true;
    
    return
end


%
max_fevals = 20000; % Maximum number of function evaluations allowed.


if ~check_params(tol_stop, 'double', true)
    disp('max_fevals is a double non-negative param')
    error_params = true;
    
    return
end


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implicit Filtering Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau= 1.0e-2; % If theta(x_k,h_k) >= -tau the direction is not a descent direction

if ~check_params(tau, 'double', true)
    disp('tau is a double non-negative param')
    error_params = true;
    
    return
end

%options = optimset('linprog','MaxIter',10000000,'TolFun',1.0e-10,'Display','none','Algorithm','dual-simplex' ); % linear programming options
options = optimset('MaxIter',10000000,'TolFun',1.0e-10,'Display','none'); % linear programming options

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimized Goldstein Line Search Options.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma_1 = 1.0e-6; %% Starting gamma F(x_k + v_k) <= F(x_k) + gamma_1 J(x_k)'v_k

if ~check_params(gamma_1, 'double', true)
    disp('gamma_1 is a double non-negative param')
    error_params = true;
    
    return
end


gamma_2 = 0; % F(x_k + alpha_k v_k) <= F(x_k) + gamma_2 alpha_k J(x_k)'v_k

if ~check_params(gamma_2, 'double', true)
    disp('gamma_2 is a double non-negative param')
    error_params = true;
    
    return
end


beta_goldstein = 0.5; % Decreasing step for Goldstein Line Search

if ~check_params(beta_goldstein, 'double', true)
    disp('beta_goldstein is a double non-negative param')
    return
end


minumum_line_search_stepsize = 1.0e-5; % Starting step for Goldstein line search

if ~check_params(minumum_line_search_stepsize, 'double', true)
    disp('minumum_line_search_stepsize is a double non-negative param')
    error_params = true;
    
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithmic Options.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialization.
%
list = 0;   % 0-4 variable: 0 if the algorithm initializes the iterate
% list with a single point; 1 if a latin hypercube sampling
% strategy is considered for initialization; 2 if
% random sampling is used; 3 if points are considered
% equally spaced in the line segment, joining the
% variable upper and lower bounds; 4 if the algorithm is
% initialized with a list provided by the optimizer;
%
user_list_size = 0;  % 0-1 variable: 1 if the user sets below the iterate
% list initial size; 0 if the iterate list initial size
% equals the problem dimension.
%
nPini = 30; % Number of points to consider in the iterate
% list initialization, when its size is defined
% by the user.
%

%
% List sorting.
%
spread_option = 1; % 0-2 variable: 1 if for each point in the current approximation
% to the Pareto frontier, each component of the objective
% function is projected in the corresponding dimension and
% points are ordered according to the largest gap between
% consecutive points, just before polling; 2 if points are
% ordered according to the largest gap between consecutive
% points in the current approximation to the Pareto frontier,
% measured using Euclidean distance. Each point has only one
% consecutive point lying in the approximation to the Pareto
% frontier (the closest one); 0 if no ordering strategy is
% considered.
%
% Directions and step size.
h_ini  = 1;    % Initial step size.
if ~check_params(h_ini, 'double', true)
    disp('h_ini is a double non-negative param')
    error_params = true;
    
    return
end


delta = 0.5;   % Decreasing factor after an unsuccessful iteration

if ~check_params(delta, 'double', true)
    disp('delta is a double non-negative param')
    error_params = true;
    return
end
