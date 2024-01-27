function [Plist,Flist,StepList,func_eval] = moif(Pareto_front,func_F,x_ini, lbound, ubound)
%
% Purpose:
%
% Function moif applies the implicit filtering method to the derivative-free
% multiobjective optimization problems with box constraints:
%
%    min F(x) = (f_1(x),f_2(x),...,f_m(x))  s.t.  lbound <=  x <= ubound
%
% where x is a real vector of dimension n. The derivatives of the functions
% f_i, i = 1,..., m, are not used.
%
% The user must provide: func_F (for F function values),
%
%
%
% Input:
%
%         Pareto_front (0-1 variable: 1 if the complete Pareto front
%                      should be computed; 0 if only one point in
%                      the Pareto front is required.)
%
%         func_F (Name of the file defining the objective function.)
%
%
%         x_ini (Initial point to be considered; only required when a file
%               is nor used for initialization, list is set equal to 0 and
%               bounds are not provided.)
%
%         lbound (Lower bounds on the problem variables. Also used for the
%                iterate list initialization.)
%
%         ubound (Upper bounds on the problem variables. Also used for the
%                iterate list initialization.)
%
%
% Output:
%
%         Plist (List of current nondominated points.)
%
%         Flist (List of function values corresponding to current
%               nondominated points.)
%
%         StepList (List of step size parameters.)
%
%         func_eval (Total number of function evaluations.)
%
%
% Functions called: func_F, parameters_moif, paretodominance (provided by the optimizer).
%
% DMS Version 0.2.
%
% Copyright (C) 2011 A. L. Custodio, J. F. A. Madeira, A. I. F. Vaz,
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
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
parameters_moif;
if error_params
    Plist = [];
    Flist = [];
    StepList = [];
    func_eval = 0;
    return
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization Step.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if output
    print_params;
    
end
%
% Define the problem size.
%
if ~isempty(lbound)
    n = size(lbound,1);
else
    if ~isempty(ubound)
        n = size(ubound,1);
    else
        if ~isempty(x_ini)
            n = size(x_ini,1);
        else
            fprintf('Error: An initial point or variable bounds should be provided.\n\n');
            return
        end
    end
end


%
% Define the initial iterate list.
%
if (isempty(lbound) || isempty(ubound) || ((sum(isfinite(lbound))+...
        sum(isfinite(ubound)))~= (2*size(lbound,1)))) && isempty(x_ini)
    fprintf('Error: An initial point or complete variable bounds should be provided. \n\n');
    return
else
    if (list == 0) || (list == 2)
        if ~isempty(x_ini)
            Pini = [x_ini];
        else
            Pini = [(lbound + ubound)/2];
        end
    else
        if ~Pareto_front
            fprintf('You cannot instantiate a list of points for a single Pareto point optimization. \n\n');
            Plist = [];
            Flist = [];
            StepList = [];
            func_eval = 0;
            return
        else
            Pini = repmat(lbound_gen,1,nPini)+repmat([0:(nPini-1)]/(nPini-1),n,1)...
                .*repmat((ubound_gen-lbound_gen),1,nPini);
        end
    end
end
%
Flist     = [];
Plist     = [];
StepList  = [];
func_eval = 0;
%
% Evaluate the initial iterate list.
%
for i=1:size(Pini,2)
    x_ini = Pini(:,i);
    %
    %  Check feasibility, without considering box projections.
    %
    feasible = 1;
    
    if ~isempty(find(ubound(:) < x_ini(:),1)) || ~isempty(find(x_ini(:) < lbound(:),1))
        feasible = 0;
    end
    
    %
    if feasible
        %
        %      Evaluate the point and store the corresponding values.
        %
        Ftemp     = feval(func_F,x_ini);
        func_eval = func_eval+1;
        
        
        if isempty(Flist)
            Flist = [Flist,Ftemp];
            Plist = [Plist,x_ini];
            StepList = [StepList,h_ini];
        else
            [pdom,index_ndom] = paretodominance(Ftemp,Flist);
            if (pdom == 0)
                Plist = [Plist(:,index_ndom),x_ini];
                Flist = [Flist(:,index_ndom),Ftemp];
                StepList = [StepList(index_ndom),h_ini];
            end
        end
    end
end
%
% Check if the iterate list is not empty.
%
if isempty(Flist)
    fprintf('Error: The optimizer did not generate a feasible point\n');
    fprintf('or the initial point provided is not feasible.\n\n');
    return
end

%
% Define the number of objective functions
%
n_obj = size(Flist,1);
%
%
halt     = 0;
iter     = 0;
iter_suc = 0;
%
% Print the iteration report header.
%
if output
    fprintf('Iteration Report: \n\n');
    if Pareto_front
        fprintf('| iter  | success | list size |       min h     |       max h     | fevals  |\n');
    else
        fprintf(['| iter  | success |',repmat(' ',1,n_obj*5), 'F', repmat(' ',1,n_obj*5), '|        h        | fevals |\n']);
    end
    print_iteration(Pareto_front, iter, Flist, StepList, 0, 0);
end
%

while (~halt)
    move      = 0;
    if (StepList(1) >= tol_stop)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     Poll Step.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %         poll = true;
        success = false;
        %         while poll
        %
        %        Generate the positive basis.
        %
        D = [eye(n) -eye(n)];
        
        nd = size(D,2);
        %
        %        Poll using the positive basis.
        %
        count_d            = 1;
        % Boolean variable which means that the current point is the
        % center of the stencil
        index_poll_center  = 1;
        
        added              = zeros(1,size(Plist,2));
        if (spread_option ~= 0)
            [Plist,Flist,StepList] = sort_gamma(Plist,Flist,StepList,true,...
                tol_stop,spread_option);
        end
        
        %
        %% Declare the stencil S(x,h) = [F_1(x+h*e_1) F_2(x+h*e_1);F_1(x+h*e_2) F_2(x+h*e_2);F_1(x-h*e_1) F_2(x-h*e_1);F_1(x-h*e_2) F_2(x-h*e_2)]
        %
        %
        %
        F_stencil = zeros(nd, n_obj);
        
        if box_projections
            P_stencil = zeros(nd,1);
        else
            P_stencil = [];
        end
        %%
        while (count_d <= nd)
            xtemp = Plist(:,1) + StepList(1) * D(:,count_d);
            %
            %           Check feasibility.
            %
            feasible = 1;
            
            %% Projection
            if box_projections
                ix_ubound = find( xtemp > ubound);
                
                if ~isempty(ix_ubound)
                    xtemp(ix_ubound) = ubound(ix_ubound);
                end
                
                ix_lbound = find( xtemp < lbound);
                
                if ~isempty(ix_lbound)
                    xtemp(ix_lbound) = lbound(ix_lbound);
                end
            end
            
            if ~isempty(find(ubound < xtemp,1)) || ~isempty(find( xtemp < lbound,1))
                feasible = 0;
            end
            
            %
            if feasible
                Ftemp = feval(func_F, xtemp);
                func_eval = func_eval + 1;
                if box_projections
                    if count_d > n
                        P_stencil(count_d) = xtemp(count_d-n);
                    else
                        P_stencil(count_d) = xtemp(count_d);
                    end
                end
                
                % Place x+h_ei into the stencil matrix
                F_stencil(count_d,:) = Ftemp';
                
                [pdom,index_ndom] = paretodominance(Ftemp,Flist);
                if (pdom == 0)
                    code_add = 1;
                    % If we are building the front, the non dominance
                    % is good :)
                    if (Pareto_front == 1)
                        success = 1;
                    end
                    
                    if index_ndom(1) == 0
                        % If we are looking for a single point, then
                        % it will be a good thing if the current point
                        % is dominated
                        if (Pareto_front == 0)
                            success  = 1;
                            code_add = 2;
                        end
                        index_poll_center = 0;
                        index_ndom(1)     = 1;
                    end
                    
                    Plist = [Plist(:,index_ndom),xtemp];
                    Flist = [Flist(:,index_ndom),Ftemp];
                    StepList  = [StepList(index_ndom),StepList(1)];
                    added = [added(index_ndom),code_add];
                end
                
            else
                % Place infinity values in the matrix
                F_stencil(count_d,:) = Inf;
            end
            count_d = count_d + 1;
        end
        %             poll = 0;
        %         end
        %
        %     Update the counter for successful iterations.
        %
        % If the coordinate-search was unsuccessfull then we try the line
        % search step
        if success
            iter_suc = iter_suc + 1;
        else
            %% Line Search Step
            if StepList(1) <= purityPercentage * max(StepList)
                % Compute h-approximated Jacobian
                J_k = compute_h_gradient(Plist(:,1), Flist(:,1), F_stencil, StepList(:,1), box_projections, P_stencil);
                
                %% LP problem
                % LP objective function
                c = [zeros(1,n) 1]; % min 0'* x + x_slack
                % constraints
                % J_k*x - J_k*x_k <= x_slack
                A = [J_k -ones(n_obj,1)];
                b = J_k * Plist(:,1);
                
                [x_hat,theta]=linprog(c,A,b,[],[],[lbound; -Inf],[ubound; +Inf],[],options);
                
                %% Take only x variable
                x_hat = x_hat(1:n);
                
                %% Compute direction
                v = x_hat-Plist(:,1);
                
                if theta >= -tau
                    success = 0;
                else
                    [alpha_k, fevals_ls, F_ls, fail_ls]=goldstein(J_k,Plist(:,1), v, func_F, Flist(:,1),minumum_line_search_stepsize, gamma_1, gamma_2, beta_goldstein);
                    
                    func_eval = func_eval + fevals_ls;
                    
                    if ~fail_ls
                        [pdom,index_ndom] = paretodominance(F_ls, Flist);
                        if (pdom == 0)
                            code_add = 1;
                            if (Pareto_front == 1)
                                success = 1;
                            end
                            if index_ndom(1) == 0
                                if (Pareto_front == 0)
                                    success  = 1;
                                    code_add = 2;
                                end
                                index_poll_center = 0;
                                index_ndom(1)     = 1;
                            end
                            Plist = [Plist(:,index_ndom), Plist(:,1) + alpha_k * v];
                            Flist = [Flist(:,index_ndom), F_ls];
                            StepList  = [StepList(index_ndom), StepList(1)];
                            added = [added(index_ndom), code_add];
                        end
                    end
                    if success
                        iter_suc = iter_suc + 1;
                    end
                end
                
            end
        end
        
        %
        %     Update the iterate list.
        %
        
        % If index_poll_center is zero, then the current stencil center is
        % dominated
        if (index_poll_center == 0)
            if (Pareto_front == 0)
                % Trick for single point optimization. The current list
                % could be composed of several points. The algorithm take
                % the first which dominates the poll center.
                index = find(added == 2);
                Plist = Plist(:,index(1));
                Flist = Flist(:,index(1));
                StepList  = StepList(index(1));
                added = added(index(1));
            else
                nPlist = size(Plist,2);
                Plist  = Plist(:,2:nPlist);
                Flist  = Flist(:,2:nPlist);
                StepList   = StepList(2:nPlist);
                added  = added(2:nPlist);
            end
        else
            if (Pareto_front == 0)
                Plist = Plist(:,1);
                Flist = Flist(:,1);
                StepList  = StepList(1);
                added = 1;
            else
                % If the poll center is not dominated, then
                % the sorting type could modify this sorting.
                nPlist = size(Plist,2);
                Plist  = [Plist(:,2:nPlist),Plist(:,1)];
                Flist  = [Flist(:,2:nPlist),Flist(:,1)];
                StepList   = [StepList(2:nPlist),StepList(1)];
                added  = [added(:,2:nPlist),1];
            end
        end
        %
        %     Update the step size parameter, if necessary.
        %
        if ~success
            StepList(logical(added))= StepList(logical(added))*delta;
        end
    else
        % The current point is not feasible due to its corresponding step <
        % tol stop.
        nPlist = size(Plist,2);
        Plist  = [Plist(:,2:nPlist),Plist(:,1)];
        Flist  = [Flist(:,2:nPlist),Flist(:,1)];
        StepList   = [StepList(2:nPlist),StepList(1)];
        move   = 1;
    end
    
    if ~move
        %
        %     Check if the stopping criteria are satisfied.
        %
        if sum(StepList >= tol_stop) == 0
            halt = 1;
        end
        if func_eval >= max_fevals
            halt = 1;
        end
        iter = iter + 1;
        %
        %     Print the iteration report.
        %
        if output
            print_iteration(Pareto_front, iter, Flist, StepList, success, func_eval)
        end
        
        %
    end
    
end
%
% Print final report in screen.
%
fprintf('\n Final Report: \n\n');
if Pareto_front
    fprintf('| #iter | #isuc   | list size |       min h     |       max h     | #fevals |\n');
else
    fprintf(['| #iter |   #isuc   |',repmat(' ',1,n_obj*5), 'F', repmat(' ',1,n_obj*5) ,'|        h        | #fevals |\n'])
end
print_iteration(Pareto_front, iter, Flist, StepList, iter_suc, func_eval)
%
%
% End moif.
