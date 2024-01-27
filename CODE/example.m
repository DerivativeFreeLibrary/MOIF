clear all; clc;

lbound = zeros(3,1);
ubound = ones(3,1);
x_ini = (lbound + ubound)./2;
Pareto_front = 1;

[Plist,Flist,StepList,func_eval] = moif(Pareto_front,@func_F,x_ini, lbound, ubound);