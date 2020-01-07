%% Control of bacteria on endogenous metabolism
close all
plot_control(1, 100)
plot_control2(15, 100)

%% Appendix: System response to Ratio F and B in input
close all

load('CONTROL_T24_FXR0_N100.mat')
for it = 1:size(control.curve.v,2)
    v_plot(it) = control.curve.v{it};
end

figure();

x   = control.curve.input_change2;
y1  = [v_plot.B1];
y2  = [v_plot.F1];

area(x, [y1; y2]'); hold on
colormap winter
xlabel('Ratio F_{in} / (F_{in} + B_{in})')
ylabel('Bacterial population (co_1)')
legend('Bacteriodetes', 'Firmicutes')