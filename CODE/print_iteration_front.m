function [] = print_iteration_front(iter, n_points, min_step, max_step, success, fevals)
%print_iteration_front print the current iteration in case of pareto front optimization.

%% Params:
%%	iter: the current iteration number
%%	n_points: is the current number of point in the list
%%	min_step: is the current minium step-size in the list;
%%	max_step: is the current maximum step-size in the list;
%%	success: 1 if the current iteration was successful, 0 otherwise
%%	fevals: current number of function evaluations.

%% Output:
%%	[]

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

print_format = '| %5d |    %2d   |   %5d   | %+13.8e | %+13.8e |  %5d  |\n';
fprintf(print_format, iter, success, n_points, min_step, max_step, fevals);

end

