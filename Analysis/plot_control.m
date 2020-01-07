function [] = plot_control(ind, n)
%PLOT_CONTROL Plot a previously generated control plot / sensitivity analysis
%   Inputs:
%   ind     - 1-8  : RYGB hypothesis plots
%             11-23: Parameters 1-13
%             24   : ratio between B_in and F_in
%             25   : k_u
%   n       - steps
% Note, the control file must be in existence already!!!!

clear v_plot1 v_plot2 x1 x2 

% NOTE - if it is the # 25, n = 100 one, only the first 82 should be
% considered
if ind == 15 && n == 100
    max = 82;
else
    max = size(control.curve.v,2);
end

load(['CONTROL_T',num2str(ind+10),'_FXR0_N',num2str(n),'.mat'])
for it = 1:max
    if ~isempty(control.curve.v{it})
        v_plot1(it) = control.curve.v{it};
    end
end
x1   = control.curve.input_change2(~isnan(control.curve.input_change(1:max)));

load(['CONTROL_T',num2str(ind+10),'_FXR1_N',num2str(n),'.mat'])
for it = 1:max
    if ~isempty(control.curve.v{it})
        v_plot2(it) = control.curve.v{it};
    end
end
x2   = control.curve.input_change2(~isnan(control.curve.input_change(1:max)));

if ind > 0
xlabels = { 'B growth rate (1/min)'
            'F growth rate (1/min)'
            'F yield'
            'B yield'
            'F recycling'
            'SCFA yield'
            'SCFA uptake'
            'Deconjugation'
            'Dehydroxylation'
            'Buffer input'
            'Fraction of HCO3- exchange'
            'Inhibition factor of F'
            'Substrate preference of F'
            'F input fraction'
            'BA synthesis'
    };

    xlabel_text = xlabels{ind};
    
elseif ind < 0
    xlabels = { 'H1c'
            'H1d'
            'H1e'
            'H2a'
            'H2b'
            'H2c'
            'H4b'
            'H5b'
            'H7'
    };
    xlabel_text = xlabels{ind+10};
end

if ind < 15
    xticks = [1e-5:9e-5 1e-4:9e-4 0.001:0.009 0.01:0.09 0.1:0.9 1:9 10:90 100:900 1e3:9e3];
    xticklabels = [1e-5 1e-4 0.001 0.01 0.1 1 10 100 1e3];
    log_scale = 0;
else
    xticks = [1e-5:9e-5 1e-4:9e-4 0.001:0.009 0.01:0.09 0.1:0.9 1:9 10:90 100:900 1e3:9e3];
    xticklabels = [0.1:0.1:1 2:1:10];
    log_scale = 1;
end

h = interaction_plot(x1, x2, v_plot1, v_plot2, xlabel_text, xticks, xticklabels, log_scale);

end

