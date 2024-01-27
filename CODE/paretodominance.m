function [pdom,index_ndom] = paretodominance(F,Flist); 
%
% Purpose:
%
%    Function paretodominance checks if the vector F satisfies a Pareto
%    dominance criterion.
%
% Input:  
%
%         F (Vector to be checked.)
%
%         Flist (List storing columnwise the objective function values 
%                at nondominated points.)
%
% Output: 
%
%         pdom (0-1 variable: 1 if the point is dominated; 0 otherwise.)
%
%         index_ndom (0-1 vector: 1 if the corresponding list point is not
%                    dominated; 0 otherwise.)
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
index_ndom    = [];
[nlist,mlist] = size(Flist);
Faux          = repmat(F,1,mlist);
index         = logical(sum(logical(Faux<Flist)));
if (sum(index) < mlist)
    pdom = 1;
else
    pdom       = 0;
    index      = sum(logical(Faux<=Flist));
    index_ndom = ~logical(index == nlist);
end
%
% End of paretodominance.