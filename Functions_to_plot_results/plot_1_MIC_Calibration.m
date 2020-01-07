function [] = plot_1_MIC_Calibration()
% PLOT_1_CALIBRATION Plot inputs, outputs and residuals (via
% plot_characterization)

%% Prepare
c       = constants_MIC();
x0      = initial_state_MIC();

%% Load optimal parameter set
load('SELECTED_OPTIMUM2')
p_opt = p_opt_adapted;

%% Plot results
tic; plot_characterization(p_opt,c,x0); toc
