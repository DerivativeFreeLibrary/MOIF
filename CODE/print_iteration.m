function [] = print_iteration(pareto_front, iter, Flist, Step_list, success, fevals)
%print_iteration calls the right function for printing the current iteration
%status
%%   Params:
%%           pareto_front: 0 if the single point iteration has to be printed, 1 otherwise
%%           iter: the current iteration number
%%           Flist: is the current Function list
%%           Step_list: is the current Step list;
%%           success: 1 if the current iteration was successful, 0 otherwise
%%           fevals: current number of function evaluations.


%%   Output:
%%           []


% MOIF v.0.1
% Copyright (C) 2017 G.Cocchi, G.Liuzzi, A.Papini, M.Sciandrone:
% This software is published under GNU General Public Licence (GPLv3). 
% Do not distribute this software under different licence, but you are allowed to use, modify and redistribute
% it under the same licence.
%
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.

if pareto_front
    n_points = size(Flist,2);
    max_step = max(Step_list);
    min_step = min(Step_list);
    print_iteration_front(iter, n_points, min_step, max_step, success, fevals)
else
    print_iteration_single(iter, Flist, Step_list, success, fevals)
end
end

