function [F] = func_F(x)
%
% Purpose:
%
%    Function func_F is an user provided function which
%    computes the value of the objective function at a
%    point provided by the optimizer.
%
% Input:  
%
%         x (Point given by the optimizer.)
%
% Output: 
%
%         F (Function value at the given point.)        
%
% DMS Version 0.2.
%
% Copyright (C) 2011 A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, 
% and L. N. Vicente.
%
% http://www.mat.uc.pt/dms
%

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
%
%
% -----------------------------------------------------------------------------
%  Block to be user modified.
% -----------------------------------------------------------------------------
f(1) = sum(x.^2);
f(2) = sum((x-5).^2);
F    = f';
% -----------------------------------------------------------------------------
%  End of block to be user modified.
% -----------------------------------------------------------------------------
%
% End of func_F.