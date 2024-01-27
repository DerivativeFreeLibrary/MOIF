%%% PRINT PARAMS %%%

disp('%%%%% MOIF Params %%%%%')

fprintf('max_h: \t %20d \n', h_ini)
fprintf('min_h: \t %20d \n', tol_stop)
fprintf('max_fevals: \t %6e \n \n', max_fevals)

disp('%%%%% Linear Programming Params %%%%% ')
fprintf('Purity Coefficient \t %10d \n', purityPercentage)
fprintf('Tau \t %26d \n \n', tau)

disp('%%%%% Goldstein Params %%%%%')
fprintf('initial step size: \t %10f \n', minumum_line_search_stepsize)
fprintf('gamma_1: \t %18d \n', gamma_1)
fprintf('gamma_2: \t %18d \n', gamma_2)
fprintf('delta: \t %26d \n \n', delta)