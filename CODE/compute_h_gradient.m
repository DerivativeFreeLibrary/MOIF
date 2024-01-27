function [ J ] = compute_h_gradient( P_center, F_center, F_stencil, h, box_projection, P_stencil)
%compute_h_gradient computes the h-approximated Jacobian.

%% Params:
%%	P_center: is the current point
%%	F_center: is the current Function list
%%	stencilValues: are F(x+/-he_i) values. It is a (2*dim) X (n_obj) matrix.
%%	h: is the approximation step;
%%	lbound: is the lower bound of the problem for the projection
%%	ubound: is the upper bound of the problem for the projection
%%	stencilPoints: Proj(x+/-h_ei) if box_projection is true else []


%% Output:
%%	J: approximated Jacobian

%% Example
%           stencilValues(i,j) = F_j(x_{mod(i,dim)+c_i * h*e_{mod(i,dim)}})
%           where c_i = 1 if i <= dim,  c_i=-1 ,otherwise
%           for example if there are two objectives and x is
%           monodimensional then
%           stencilValues = [F_1(x+h) F_2(x+h); F_1(x-h) F_2(x-h)]

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



n = length(P_center);
n_obj = size(F_center,1);

% The Jacobian matrix is instantiated
J = zeros(n_obj, n);
infinity_vector = Inf*ones(1,n_obj);
for i=1:n
    % F(x_i+he_i)
    F_i_positive = F_stencil(i,:);
    
    % F(x_i-he_i)
    F_i_negative = F_stencil(i+n,:);
    
    if ~box_projection
        % Infeasibility, then finite difference will be left or right
        feasible = true;
        if sum(F_i_positive == infinity_vector) == n_obj
            F_i_positive = F_center';
            feasible = false;
        end
        if sum(F_i_negative == infinity_vector) == n_obj
            F_i_negative = F_center';
            feasible = false;
        end
        
        if feasible
            step = 2*h;
        else
            step = h;
        end
        
        
    else
        step = P_stencil(i)-P_stencil(i+n);
    end

    % In case of infeasibility with respect to left and right moves, then
    % the Jacobian will be zero
    J(:,i) = (F_i_positive-F_i_negative)./step;
end

end

