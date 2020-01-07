function [h_out] = plot_hypotheses(reg, n, h, use_colors)
tic
%% Figure prep
if nargin < 3
    h_out = figure('units', 'normalized', 'Position', [0 0 .8 .8]);
    use_colors = 1;
else
    figure(h);
    h_out = h;
end
ys              = 4;
xs              = 8;

% Hypothesis colors
g               = [0.7 0.7 0.7];
if use_colors
cs              = { [0 0.5 0]       % H1c, dark green
    [0 0.7 0]       % H1d, dark green
    [0 0.9 0]       % H1e, dark green
    [0.5 1 0]       % H2a, yellow-green
    [0.5 0.8 0]     % H2b, yellow-green
    [0.5 0.6 0]     % H2c, yellow-green
    [0.9 0.5 0]     % H4b, orange
    [0.9 0 0]       % H5b, red
    [0 0 0]         % H7,  black
    };
else
    cs              = {g       % H1c, dark green
    g      % H1d, dark green
    g      % H1e, dark green
    g       % H2a, yellow-green
    g     % H2b, yellow-green
    g     % H2c, yellow-green
    g     % H4b, orange
    g       % H5b, red
    g         % H7,  black
    };
end



%% Load simulations
input_change    = zeros(6,n);
input_change2   = zeros(6,n);
normal_x        = zeros(8,1);
normal_y        = zeros(4,1);

for it = 1:8
    clear v_plot
    
    load(['CONTROL_T',num2str(it),'_FXR',num2str(reg),'_N',num2str(n),'.mat'])
    
    for it2 = 1:n
        if ~isnan(control.curve.input_change(it2))
            v_plot(it2)     = control.curve.v{it2};
        end
    end
    
    v_pre           = control.ref.v(end);
    input_change    = control.curve.input_change;
    input_change2   = control.curve.input_change2;
    
    if it == 5 % Normal value for pH is 6.5, or a change of -0.5
        normal_x(it)    = interp1(input_change(~isnan(input_change)), input_change2(~isnan(input_change)), -0.5);
    else
        normal_x(it)    = interp1(input_change(~isnan(input_change)), input_change2(~isnan(input_change)), 1);
    end
    
    %% Plot
    ind = ~isnan(input_change);
    
    subplot(ys, xs, 0*8+it)
    y_normal(1) = v_pre(end).BA_pool;
    plot([normal_x(it) normal_x(it)], [1e3, 1e5],               	'k', 'LineWidth', 0.5); hold on
    plot([-1e3, 1e3],  [y_normal(1) y_normal(1)],   'k', 'LineWidth', 0.5); hold on
    plot(input_change2(ind),[v_plot.BA_pool],                                    '-', 'Color', cs{it}, 'MarkerFaceColor', cs{it}, 'MarkerSize', 10, 'LineWidth', 2); hold on 
%     plot(input_change(it,end),[v_plot(end).BA_pool],                                   'o', 'Color', cs{it}, 'MarkerFaceColor', cs{it}, 'MarkerSize', 10, 'LineWidth', 2); hold on
    xlim([min(input_change2) max(input_change2)]); 
    
    subplot(ys, xs, 1*8+it)
    y_normal(2) = [v_pre(end).tr_F_6_5]./[v_pre(end).tr_B_6_5];
    plot([normal_x(it) normal_x(it)], [-1e3, 1e3],               	'k', 'LineWidth', 0.5); hold on
    plot([-1e3, 1e3],  [y_normal(2) y_normal(2)],   'k', 'LineWidth', 0.5); hold on
    plot(input_change2(ind),[v_plot.tr_F_6_5]./[v_plot.tr_B_6_5],                '-', 'Color', cs{it}, 'MarkerFaceColor', cs{it}, 'MarkerSize', 10, 'LineWidth', 2); hold on 
    xlim([min(input_change2) max(input_change2)]); 
    
    subplot(ys, xs, 2*8+it)
    y_normal(3) = [v_pre(end).pH1];
    plot([normal_x(it) normal_x(it)], [-1e3, 1e3],               	'k', 'LineWidth', 0.5); hold on
    plot([-1e3, 1e3],  [y_normal(3) y_normal(3)],   'k', 'LineWidth', 0.5); hold on
    plot(input_change2(ind),[v_plot.pH1],                                        '-', 'Color', cs{it}, 'MarkerFaceColor', cs{it}, 'MarkerSize', 10, 'LineWidth', 2); hold on 
    xlim([min(input_change2) max(input_change2)]); 
    
    subplot(ys, xs, 3*8+it)
    y_normal(4) = 100*[v_pre(end).fS];
    plot([normal_x(it) normal_x(it)], [-1e3, 1e3],               	'k', 'LineWidth', 0.5); hold on
    plot([-1e3, 1e3],  [y_normal(4) y_normal(4)],   'k', 'LineWidth', 0.5); hold on
    plot(input_change2(ind),100*[v_plot.fS],                                     '-', 'Color', cs{it}, 'MarkerFaceColor', cs{it}, 'MarkerSize', 10, 'LineWidth', 2); hold on 
    xlim([min(input_change2) max(input_change2)]); 
       
end

%% Legends
% 
rev   = [1 1 1 1 0 1 1 1];
names = {'Carb input (H1c)' 'RS input (H1d)' 'Nutrient input (H1e)' 'BA cycling (H2a)' 'Buffer pH change (H2b)' 'Buffer input (H2c)' 'Colon speed (H4b)'  'Frac colon input (H5b)'};
for it = 1:length(names)
    subplot(ys,xs,it)
    title(names{it})
    
    if it == 1
        ylabel('BA pool (\mumol)')
    end
    ylim([2500 11000])
    
    if rev(it)
        set(gca, 'XDir', 'reverse')
    end
    
    subplot(ys,xs,xs*1+it)
    if it == 1
        ylabel('Fecal F/B ratio')
    end
    ylim([0 4])
    
    if rev(it)
        set(gca, 'XDir', 'reverse')
    end
    
    subplot(ys,xs,xs*2+it)
    if it == 1
        ylabel('Proximal pH')
    end
    ylim([5 7])
    
    if rev(it)
        set(gca, 'XDir', 'reverse')
    end
    
    subplot(ys,xs,xs*3+it)
    switch it 
        case 1
            ylabel('Secondary BA (%)')
            xlabel('Non-RS input (g/day)', 'FontSize', 8)
        case 2
            xlabel('Resistant starch input (g/day)', 'FontSize', 8)
        case 3
            xlabel('Carbohydrate input (g/day)', 'FontSize', 8)
        case 4
            xlabel('Bile acid cycling (pools/day)', 'FontSize', 8)
        case 5
            xlabel('Buffer pH', 'FontSize', 8)
        case 6
            xlabel('Buffer input (mol/day)', 'FontSize', 8)
        case 7
            xlabel('Colon transit (1/min)', 'FontSize', 8)
        case 8
            xlabel('Fractional colon input (%)', 'FontSize', 8)
    end
    ylim([0 100])
    
    if rev(it)
        set(gca, 'XDir', 'reverse')
    end
end


if use_colors
    if reg == 0
        suptitle('Microbiota model, no FXR regulation');
    elseif reg == 1
        suptitle('Microbiota model, FXR regulation');
    end
end


toc
