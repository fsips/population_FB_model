function [h] = interaction_plot(x1, x2, v_plot1, v_plot2, xlabel_text, xticks, xticklabels, log_scale)
%INTERACTION_PLOT Generate a control plot / sensitivity analysis
%
%   Inputs:
%
%       x1, x2, v_plot1, v_plot2 - contain the xs and ys for the left (1) and
%           right (2) control plots
%
%       xlabel_text, xticks, xticklabels - plot annotation
%
%       log_scale - if the BA pool becomes very large, choose to not linearly
%           scale BA pool and co1 bile acids
% 
% Figure handle h is returned

h = figure('Units', 'normalized', 'Position', [0.1 0 0.4 1]);

c_B     = [1/5 4/5 1];
c_F     = [2/5 4/5 7/10];
c_C1    = [2/5 4/5  3/10];
c_C2    = [3/5 9/10 5/10];
c_C3    = [4/5 1    7/10];
cm_C    = [c_C1; c_C2; c_C3];

labs = {};
for it = 1:length(xticks)
    if sum (xticks(it) == xticklabels)
        labs{it} = num2str(xticks(it));
    else
        labs{it} = '';
    end
end

for it = 1:2
    d = it-1;
    switch it
        case 1
            v_plot = v_plot1;
            x      = x1;
        case 2
            v_plot = v_plot2;
            x      = x2;
    end
    
    for it2 = 1:4
        switch it2
            case 1
                
                subplot(4,2,1+d)
                
                y1  = [v_plot.B1];
                y2  = [v_plot.F1];
                
                h = area(x, [y1; y2]'); hold on
                
                set(h(1), 'FaceColor', c_B)
                set(h(2), 'FaceColor', c_F)
                
                ylabel('Bacterial population (co_1)')
                if d
                    legend('B', 'F', 'Box', 'off', 'Color', 'none', 'Location', 'SouthEast')
                    title('FXR regulation')
                else
                    title('No FXR regulation')
                end
                
            case 2
                
                subplot(4,2,3+d)
                
                y1  = [v_plot.PU_1];
                y2  = [v_plot.PC_1];
                y3  = [v_plot.SU_1];
                y4  = [v_plot.SC_1];
                
                h = area(x, [y1;y2;y3;y4]'); hold on
                
                set(h(1),'FaceColor',[190 186 205]/255)
                set(h(2),'FaceColor',[146 137 204]/255)
                set(h(3),'FaceColor',[190 172 193]/255)
                set(h(4),'FaceColor',[175 113 188]/255)
                
                if log_scale 
                else
                    ylim([0 700])
                end
                ylabel('Bile acids in co_1')
                if d
                    legend('PU', 'PC', 'SU', 'SC', 'Box', 'off', 'Color', 'none', 'Location', 'SouthEast')
                end
                
            case 3
                
                subplot(4,2,5+d)
                
                y1  = [v_plot.pH1];
                
                plot(x, y1,'-', 'LineWidth', 2, 'Color', [85 88 165]/255 ); hold on
                
                ylim([5.5 6.5])
                ylabel('pH in co_1')
                
            case 4
                
                subplot(4,2,7+d)
                
                y1  = [v_plot.BA_pool];
                
                plot(x, y1,'-', 'LineWidth', 2, 'Color', [190 82 131]/255 ); hold on

                if log_scale 
                    set(gca, 'YScale', 'log')
                else
                    ylim([0 1e4])
                end
                ylabel('BA pool (\mu mol)')
                
                xlabel(xlabel_text)
                
        end
        
        
        xlim([min(x) max(x)])
        set(gca, 'XScale', 'log')
        set(gca, 'XTick', xticks)
        set(gca, 'XTickLabel', labs)
    end
end