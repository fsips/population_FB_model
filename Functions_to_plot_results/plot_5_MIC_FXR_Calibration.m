function [] = plot_5_MIC_FXR_Calibration()
% PLOT_5_MIC_FXR_CALIBRATION Plots an overview of the FXR calibration that
% yielded a parameter for the FXR feedback gain

%% LOAD Calibration
load('RES_Calibration.mat')

%% Plot calibration
h = figure();
set(h, 'Position', [595   558   865   342]);

% First: Surface plot of synthesis change vs disappearance parameter and
% FXR gain parameter
subplot(1,2,1)
surf(log10(synth./c.ku))

ylabel('FXR regulation')
xlabel('BA disappearance')
zlabel('Fold change of BA synthesis')

n = 25;
ind = 1:4:n;
set(gca, 'YTick', ind)
set(gca, 'YTickLabel', fxr(ind))
set(gca, 'XTick', ind)
set(gca, 'XTickLabel', dis(ind))
set(gca, 'ZTick', 0:5)
set(gca, 'ZTickLabel', 10.^[0:5])

% Then: The results with the maximal value of the disappearance, and the
% point of the parameter l_FXR
subplot(1,2,2)
loglog(fxr, synth(:,end)./c.ku, 'Color', [0 0.4 0.8], 'LineWidth', 2); hold on
plot(l_FXR, 10, 'sk', 'MarkerFaceColor', [0 0 0])

xlabel('FXR regulation')
ylabel('Fold change of BA synthesis')
xlim([10^-6 10^2])