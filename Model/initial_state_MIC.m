function [x0] = initial_state_MIC()

x0          = 1e-12*ones(62,1); 
x0(1:5)     = 0.2;
x0(6)       = 1;
x0(11)      = 1;