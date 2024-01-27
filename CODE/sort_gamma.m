function [Psort,Fsort,alfasort] = sort_gamma(P,F,alfa,stop_alfa,tol_stop,spread_option);
%
% Purpose:
%
%    Function sort_gamma sorts a list of points according to the largest
%    gap between two consecutive points in the list considered.
%
% Input:
%
%         P (List of points.)
%
%         F (Function values corresponding to the points in the list.)
%
%         alfa (Step size parameters corresponding to the points in the list.)
%
%         stop_alfa (0-1 variable, which indicates if there is a stopping criterion
%                    based on the step size.)
%
%         tol_stop (Lowest value allowed for the step size parameter.)
%
%         spread_option (0-2 variable, which indicates the type of ordering
%                      strategy to consider)
%
% Output:
%
%         Psort (Sorted list of points.)
%
%         Fsort (Function values corresponding to the points in the sorted list.)
%
%         alfa (Step size parameters corresponding to the points in the
%         sorted list.)
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

Psort    = P;
Fsort    = F;
alfasort = alfa;
if (spread_option ~= 0)
    if stop_alfa
        index1 = find(alfa>tol_stop);
    else
        index1 = [1:size(P,2)];
    end
    Paux   = P(:,index1);
    Faux   = F(:,index1);
    alfaux = alfa(index1);
    nPaux  = size(Paux,2);
    if (nPaux >= 2)
        [alfaux,index] = sort(alfaux,'descend');
        Paux           = Paux(:,index);
        Faux           = Faux(:,index);
        if (spread_option == 1)
            mPaux = size(Faux,1);
            Mdist = -inf*ones(2*mPaux,nPaux);
            for j = 1: mPaux
                [FProj,index] = sort(Faux(j,:));
                for i = 1:nPaux
                    if (i == 1)
                        Mdist(2*j-1,index(1)) = FProj(2)-FProj(1);
                    else
                        if (i == nPaux)
                            Mdist(2*j-1,index(nPaux)) = FProj(nPaux)-FProj(nPaux-1);
                        else
                            Mdist(2*j-1,index(i)) = FProj(i)-FProj(i-1);
                            Mdist(2*j,index(i))   = FProj(i+1)-FProj(i);
                        end
                    end
                end
            end
            [aux,index] = sort(max(Mdist),'descend');
        else
            Mdist = inf*ones(nPaux,nPaux);
            i     = 1;
            while (i<= nPaux)
                for j = (i+1):nPaux
                    Mdist(i,j) = norm(Faux(:,i)-Faux(:,j),2);
                    Mdist(j,i) = Mdist(i,j);
                end
                i = i+1;
            end
            [aux,index] = sort(min(Mdist),'descend');
        end
        Paux        = Paux(:,index);
        Faux        = Faux(:,index);
        alfaux      = alfaux(index);
        if stop_alfa
            index1 = find(alfa<=tol_stop);
        else
            index1 = [];
        end
        Psort    = [Paux,P(:,index1)];
        Fsort    = [Faux,F(:,index1)];
        alfasort = [alfaux,alfa(index1)];
    else
        if stop_alfa
            index1 = find(alfa<=tol_stop);
        else
            index1 = [];
        end
        Psort    = [Paux,P(:,index1)];
        Fsort    = [Faux,F(:,index1)];
        alfasort = [alfaux,alfa(index1)];
        
    end
end
%
% End of sort_gamma.