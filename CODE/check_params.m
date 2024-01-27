function [ is_correct ] = check_params( param, type, non_negative)
%check_params checks if the current param type and sign are correct. 

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

is_correct = true;
current_type = class(param);

if ~strcmp(type, current_type)
    is_correct= false;
else
    if strcmp(type,current_type) & non_negative
        if param < 0
            is_correct = false;
        end
    end
end

end

