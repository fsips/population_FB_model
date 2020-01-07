function [] = plot_2_MIC_Dynamics()

%% Prepare
c0      = constants_MIC();
x0      = initial_state_MIC();
load('SELECTED_OPTIMUM2');
p_opt = p_opt_adapted;
use_fxr = 1;

global bt
bt = tic;
options = odeset('NonNegative',1, 'AbsTol', 1e-9, 'RelTol', 1e-6);
[t0,s0] = ode15s(@ode_MIC,[0 1e5*24*60],x0,options,p_opt,c0,[]);

for it = 1:size(s0,1)
    v0(it) = fluxes_MIC(t0(it), s0(it,:),p_opt,c0);
end

%% Implement FXR
if use_fxr
    load('RES_Calibration.mat')
    c0.FXR   = l_FXR;
    c0.cyc0  = v0(end).cyc;
end

%% Dietary change
c1       = c0;
c1.RSin  = c0.RSin*3;

global bt
bt = tic;
[t1,s1] = ode15s(@ode_MIC,t0(end)+[0 1e5*24*60],s0(end,:),options,p_opt,c1,[]);

for it = 1:size(s1,1)
    v1(it) = fluxes_MIC(t1(it), s1(it,:),p_opt,c1);
end

%% Bowel emptying
c2       = c0;
x02      = 0.05*s0(end,:);
x02(56)  = s0(end, 56);

global bt
bt = tic;
[t2,s2] = ode15s(@ode_MIC,t0(end)+[0 1e5*24*60],x02,options,p_opt,c2,[]);

for it = 1:size(s2,1)
    v2(it) = fluxes_MIC(t2(it), s2(it,:),p_opt,c2);
end

%% Antibiotics
c3a       = c0;
x03a      = s0(end,:);
x03a(6:10)  = 0.01*s0(end,6:10);

global bt
bt = tic;
[t3a,s3a] = ode15s(@ode_MIC,t0(end)+[0 1e5*24*60],x03a,options,p_opt,c3a,[]);

for it = 1:size(s3a,1)
    v3a(it) = fluxes_MIC(t3a(it), s3a(it,:),p_opt,c3a);
end

% c3b       = c0;
% x03b      = s0(end,:);
% x03b(11:15)  = 0.01*s0(end,11:16);
% 
% global bt
% bt = tic;
% [t3b,s3b] = ode15s(@ode_MIC,t0(end)+[0 1e5*24*60],x03b,options,p_opt{1},c3b,[]);
% 
% for it = 1:size(s3b,1)
%     v3b(it) = fluxes_MIC(t3b(it), s3b(it,:),p_opt{1},c3b);
% end


%% Plot

figure()

for it = 1:3
    switch it
        case 1
            t_dyn = t1;
            v_dyn = v1;
            dp = 0;
        case 2
            t_dyn = t2;
            v_dyn = v2;
            dp = 3;
        case 3
            t_dyn = t3a;
            v_dyn = v3a;
            dp = 6;
    end
    
    subplot(3,3,1+dp)
    plot([t0;t_dyn]/(24*60) - t_dyn(1)/(24*60), [[v0.F1]';[v_dyn.F1]'], 'LineWidth', 2, 'Color', [2/5 4/5 7/10]); hold on
    plot([t0;t_dyn]/(24*60) - t_dyn(1)/(24*60), [[v0.B1]';[v_dyn.B1]'], 'LineWidth', 2, 'Color', [1/5 4/5 1]); hold on
    xlim([-1 5])
    ylim([0 15]*((2.5e11 * (100/(60*24)) * 0.9)/1e12 / 0.25))
    if it == 3
        xlabel('Time (days)');
    end
    ylabel('Bacterial content of co_1 (10^{12} #)')
    
    subplot(3,3,2+dp)
    semilogy([t0;t_dyn]/(24*60) - t_dyn(1)/(24*60), [[v0.RS_1]';  [v_dyn.RS_1]'],   'LineWidth', 2, 'Color', [2/5 4/5  3/10]); hold on
    plot([t0;t_dyn]/(24*60)     - t_dyn(1)/(24*60), [[v0.C_1]';   [v_dyn.C_1]'],    'LineWidth', 2, 'Color', [3/5 9/10 5/10]); hold on
    plot([t0;t_dyn]/(24*60)     - t_dyn(1)/(24*60), [[v0.SCFA_1]';[v_dyn.SCFA_1]'], 'LineWidth', 2, 'Color', [4/5 1    7/10]); hold on
    xlim([-1 5])
    ylim([10^-8 10^0])
    if it == 3
        xlabel('Time (days)');
    end
    ylabel('Carbohydrate content of co_1 (mol)')
    
    subplot(3,3,3+dp)
    plot([t0;t_dyn]/(24*60) - t_dyn(1)/(24*60), [[v0.PU_1]';[v_dyn.PU_1]'], 'LineWidth', 2, 'Color', [190 186 205]/255); hold on
    plot([t0;t_dyn]/(24*60) - t_dyn(1)/(24*60), [[v0.PC_1]';[v_dyn.PC_1]'], 'LineWidth', 2, 'Color', [146 137 204]/255); hold on
    plot([t0;t_dyn]/(24*60) - t_dyn(1)/(24*60), [[v0.SU_1]';[v_dyn.SU_1]'], 'LineWidth', 2, 'Color', [190 172 193]/255); hold on
    plot([t0;t_dyn]/(24*60) - t_dyn(1)/(24*60), [[v0.SC_1]';[v_dyn.SC_1]'], 'LineWidth', 2, 'Color', [175 113 188]/255); hold on
    xlim([-1 5])
    ylim([0 400])
    if it == 3
        xlabel('Time (days)');
    end
    ylabel('Bile acid content of co_1 (\mumol)')
end