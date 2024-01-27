function [ alpha_k,fevals, f_new, fail] = goldstein(J_k, x_k, v_k, func_F, f_xk, start_alpha, gamma_1, gamma_2,  beta)
%goldstein applies a Goldstien-type line search thorugh a (hopefully)
%descent direction

%% Params
%%	J_k: approximation of Jacobian computed in x_k
%%	x_k: current point
%%	v_k: approximated descent direction in x_k
%%	func_F: name of the function who computes objective function values
%%	lb: lower bound vector
%%	ub: upper bound vector
%%	start_alpha: line search minimum stepsize
%%	gamma_1: initial step approval
%%	gamma_2: following step approval
%%	delta: alpha_k = alpha_k/delta

%% Output
%%	alpha_k: stepsize found by algorithm
%%	fevals: number of function evaluations done by algorithm
%%	f_new: f(x_k+alpha_k*v_k)
%%	fail: boolean value which represent a line search failure

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
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

n_obj = size(f_xk,1);
fail = 0;

alpha = start_alpha;

f_curr = feval(func_F, x_k + alpha * v_k);

fevals = 1;

goldsteinCond = f_curr - (f_xk + gamma_1 * alpha * J_k * v_k);

%% If the descent direction is not verified the algorithm stops
if  length(find(goldsteinCond<=0))<n_obj
    fail = true;
    alpha_k = alpha;
    f_new = f_xk;
    return;
end

%% We set the actual best function value as f(xk + start_alpha * v_k)
f_best = f_curr;
%% If the descent direction is verified then the next point is obtained with a greater step
alpha = alpha/beta;
f_curr = feval(func_F, x_k + alpha * v_k);

%% Trying to extend the step in order to minimize the objectives

goldsteinCond = f_curr - ( f_best + gamma_2*alpha*J_k*v_k);

while length(find(goldsteinCond<0)) == n_obj && alpha <= 1
    alpha = alpha/beta;
    %% If alpha is greater than 1 (we are outside the box) we break the loop
    if alpha > 1
        break
    end
    %% setting the f_best value as f_curr
    f_best = f_curr;
    %% compute the next f_curr
    f_curr = feval(func_F,x_k+alpha*v_k);
    fevals = fevals + 1;
    
    goldsteinCond = f_curr - (f_best + gamma_2 * alpha * J_k * v_k);

end

%% If the above loop finished because of alpha too high, we project the point on the box

if alpha > 1
    f_curr = feval(func_F, x_k + v_k);
    fevals = fevals + 1;
    goldsteinCond = f_curr-(f_best+gamma_2*J_k*v_k);
    if length(find(goldsteinCond<0)) == n_obj
        alpha_k = 1;
        f_new = f_curr;
        return
    end
end

f_new = f_best;
alpha_k = alpha * beta;


end

