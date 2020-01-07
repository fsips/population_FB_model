function [] = generate_control_plot(use_fxr, type, n)
%GENERATE_CONTROL_PLOT2 Generate a control plot / sensitivity analysis
%   Inputs:
%   use_fxr - 0    : model wihtout FXR feedback
%             1    : model with FXR feedback
%   type    - 1-8  : RYGB hypothesis plots
%             11-23: Parameters 1-13
%             24   : ratio between B_in and F_in
%             25   : k_u
%   n       - Number of steps in the control plot. Note that below, the
%             the range sometimes depends on the number of steps.\
% 
% RESULTS ARE SAVED, via: save(['CONTROL_T', num2str(type), '_FXR', num2str(use_fxr), '_N', num2str(n)], 'control');

disp(['CONTROL_T', num2str(type), '_FXR', num2str(use_fxr), '_N', num2str(n)]);

%% Prepare
c0      = constants_MIC();
x0      = initial_state_MIC();
days    = 1e5;
options = odeset('NonNegative',1, 'AbsTol', 1e-9, 'RelTol', 1e-6);

%% Load optimization results
load('SELECTED_OPTIMUM2.mat');
p0      = p_opt_adapted;

%% Initial simulation
global bt;
bt = tic;
[t_pre,s_pre]   = ode15s(@ode_MIC,[0 days*24*60],x0,options,p0,c0);
v_pre           = fluxes_MIC(t_pre(size(s_pre,1)), s_pre(size(s_pre,1),:),p0,c0);

control.ref.v = v_pre;
control.ref.s = s_pre;
control.ref.t = t_pre;
control.ref.x0 = x0;
control.ref.p = p0;
control.ref.c = c0;

disp('Control simulation completed.');

%% Implement FXR
if use_fxr
    load('RES_Calibration.mat')
    c0.FXR   = l_FXR;
    c0.cyc0  = v_pre(end).cyc;
end

%% Control plot
for it = 1:n
    %% Repair
    c = c0;
    p = p0;
    
    switch type
        
        % Hypotheses
        case 1  % 1st, H1c, dark green (Carbohydrate input)
            c.Cin                   = c.Cin./(1+(it-n./2)./n*1.75);
            input_change(it)      = 1./(1+(it-n./2)./n*1.75);
            input_change2(it)     = c0.Cin./(1+(it-n./2)./n*1.75).*24*60*162.1;
            
        case 2  % 2nd, H1d, dark green (RS input)
            c.RSin                  = c.RSin/(1+(it-n./2)./n*1.75);
            input_change(it)      = 1/(1+(it-n./2)./n*1.75);
            input_change2(it)     = c0.RSin/(1+(it-n./2)./n*1.75)*24*60*162.1;
            
        case 3  % 3rd, H1e, dark green
            c.Cin                   = c.Cin/(1+(it-n./2)./n*1.75);
            c.RSin                  = c.RSin/(1+(it-n./2)./n*1.75);
            input_change(it)      = 1/(1+(it-n./2)./n*1.75);
            input_change2(it)     = c0.Cin/(1+(it-n./2)./n*1.75) *24*60*162.1 +c0.RSin/(1+(it-n./2)./n*1.75) *24*60*162.1;
            
        case 4  % 4th, H2a, BA cycling
            c.rel_bO                = c.rel_bO/(1+(it-n/2)/n*1.3);
            input_change(it)      = 1/(1+(it-n/2)/n*1.3);
            input_change2(it)     = c0.rel_bO/(1+(it-n/2)/n*1.3)*24*60;
            
        case 5  % 5th, H2b, raised pH
            c.raise_pH              = (it-n*2/3)/(n/4.5);
            input_change(it)      = (it-n*2/3)/(n/4.5);
            input_change2(it)     = 7+(it-n*2/3)/(n/4.5);
            
        case 6  % 6th, H2c, buffer input
            p(10)                    = p0(10)/(1+(it-n/2)/n*1.3);
            input_change(it)      = 1/(1+(it-n/2)/n*1.3);
            input_change2(it)     =p0(10)/(1+(it-n/2)/n*1.3)*24*60;
            
        case 7  % 7th, H4b, slower colonic transit
            c.k_co                  = (1-(it-n/2)/n/2)*c.k_co;
            input_change(it)      = (1-(it-n/2)/n/2);
            input_change2(it)     =(1-(it-n/2)/n/2)*c0.k_co;
            
        case 8   % 8th, H5b, fractional colonic input
            c.rel_cI                = c.rel_cI/(1+(it-n/2)/n*1.3);
            input_change(it)      = 1/(1+(it-n/2)/n*1.3);
            input_change2(it)     = c0.rel_cI/(1+(it-n/2)/n*1.3)*100;
            
            % Control plots for parameters
        case {11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23}
            it_p = type-10;
            p(it_p)               = p0(it_p)*10.^(-3+6*((it-1)/(n-1)));
            input_change(it)      = 10.^(-3+6*((it-1)/(n-1)));
            input_change2(it)     = p0(it_p)*10.^(-3+6*((it-1)/(n-1)));
            
            % Control plots for bacterial input
        case 24
            bact_in               = c0.Fin + c0.Bin;
            c.Fin                 = it / (n+1) * bact_in;
            c.Bin                 = (n+1-it)/ (n+1) * bact_in;
            input_change(it)      = it / (n+1);
            input_change2(it)      = it / (n+1);
            
            % Control plots for bile acid synthesis parameter ku
        case 25
            if n < 20
                c.ku                 = c0.ku*2.^(-1+2*((it-1)/(n-1)));
                input_change(it)      = 2.^(-1+2*((it-1)/(n-1)));
                input_change2(it)     = c0.ku*2.^(-1+2*((it-1)/(n-1)));
            else
                c.ku                 = c0.ku*2.^(-3+6*((it-1)/(n-1)));
                input_change(it)      = 2.^(-3+6*((it-1)/(n-1)));
                input_change2(it)     = c0.ku*2.^(-3+6*((it-1)/(n-1)));
            end
            
        otherwise
            disp('Control plot type unknown')
            
    end
    %% Simulation
    if type == 25
        if it<83 %% Dirty fix for the fact that # 25 , fxr0 gets stuck at 83....
        bt              = tic;
            [t,s]           = ode15s(@ode_MIC,t_pre(end)+[0 days*24*60],s_pre(end,:),[],p,c);
        end
    else
        bt              = tic;
        [t,s]           = ode15s(@ode_MIC,t_pre(end)+[0 days*24*60],s_pre(end,:),[],p,c);
    end
    disp(['Iteration ', num2str(it), ' of ',num2str(n), ' completed.']);
    
    if t(end) == t_pre(end)+days*24*60
        v      = fluxes_MIC(t(size(s,1)),s(size(s,1),:),p,c);
        control.curve.v {it} = v;
    else
        input_change(it) = NaN;
    end
    
    %% Admin
    control.curve.input_change = input_change;
    control.curve.input_change2 = input_change2;
    
end

%% Save results
save(['CONTROL_T', num2str(type), '_FXR', num2str(use_fxr), '_N', num2str(n)], 'control');