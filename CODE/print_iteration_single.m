function [] = print_iteration_single(iter, f_curr, h, success, fevals)
%print_iteration_single print the current iteration in case of single point optimization.
%% Params:
%%	iter: the current iteration number
%%	fcurr: is the current function value
%%	h: is the current step;
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


n_obj = length(f_curr);
print_format = strcat('| %5d |    %4d |', repmat(' %9f ',1,n_obj),'| %+13.8e | %5d  |\n');
fprintf(print_format, iter, success, f_curr', h, fevals);

end

