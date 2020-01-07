
%% Preparation

% Define constants, parameters, initial states
c       = constants_MIC();

load('SELECTED_OPTIMUM2.mat')
p       = p_opt_adapted;

x0      = initial_state_MIC();
days    = 1e3;
options = odeset('NonNegative',1, 'AbsTol', 1e-9, 'RelTol', 1e-6);

% Run an initial simulation, to determine values for FXR regulation
global bt;
bt = tic;
[t_pre,s_pre] = ode15s(@ode_MIC,[0 days*24*60],x0,options,p,c);
for it = 1:size(s_pre,1)
    v_pre(it)  = fluxes_MIC(t_pre(it), s_pre(it,:),p,c);
end

c.cyc0  = v_pre(end).cyc;

%% Run for various values of c.FXR
n       = 25;
fxr     = 10.^[-6:(2+6)/(n-1):2];
dis     = 10.^[-6:(2+6)/(n-1):2];
days    = 50;

for it1 = 1:length(fxr)
    for it2 = 1:length(dis)
        
        (it1-1)*length(fxr)+it2
        
        % And fxr regulation to the required amount
        c.FXR           = fxr(it1);
        
        % Set bile acid uptake and cI so that bile acids are removed
        c.dis           = dis(it2);
        
        bt              = tic;
        [t,s]           = ode15s(@ode_MIC,t_pre(end)+[0 days*24*60],s_pre(end,:),options,p,c);
        
        if t(end) == t_pre(end)+days*24*60
            v_SS            = fluxes_MIC(t(end),s(end,:),p,c);
            synth(it1,it2)  = v_SS.ku;
            
        else
            synth(it1,it2) = NaN;
        end
    end
end

l_FXR = interp1(synth(:,end)./c.ku, fxr, 10);

save('RES_Calibration.mat', 'synth', 'fxr', 'dis', 'c', 'l_FXR')

