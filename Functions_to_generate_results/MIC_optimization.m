function [] = MIC_optimization(r_set)
% Run 1: Optimization

%% Preparation
c       = constants_MIC();
load('p0.mat')
p0      = p_opt{1};
x0      = initial_state_MIC();
% plot_characterization(p0,c,x0);

%% Optimization options
options = optimset('Display', 'Iter', 'TolFun', 1e-12, 'TolX', 1e-12, 'TypicalX', p0);
lb      = ones(13,1)*-Inf;
hb      = ones(13,1)*Inf;
hb(1)   = -1;
hb(2)   = -1;
hb(6)   = 0;
hb(11)  = 0;
hb(13)  = 0;
n       = 5;
n2      = 25;
sc      = 1;


%% Optimization

% Set options, low and high bounds
E_best  = Inf;

% Set random setting
rng(r_set)
for it = 1:n
    disp(['Optimization ', num2str(it), ' of ', num2str(n)]);
    
    for it2 = 1:n2
        samples{it,it2} = p0+p0.*(randn(size(p0))/sc);
        samples_E{it,it2} = sum((OF_MIC(samples{it,it2},c,x0)).^2);
    end
    
    curr_E = [samples_E{it,:}];
    min_E  = find(curr_E == min(curr_E));
    
    [p_opts{it} E_opts(it)] = lsqnonlin(@(p_curr)OF_MIC(p_curr,c,x0), samples{it,min_E}, 10.^lb, 10.^hb, options);
    
    if E_opts(it) < E_best
        p_curr = p_opts{it};
        E_curr = E_opts(it);
        save(['CURR_BEST_',num2str(r_set),'.mat'], 'p_curr', 'E_curr');
        E_best = E_opts(it);
    end  
end

loc     = E_opts == min(E_opts);
p_opt   = p_opts(loc);
E_opt   = E_opts(loc);

% Save results
save(['RES_Optimization_n', num2str(n),'_r',num2str(r_set),'.mat'], 'p_opt', 'E_opt', 'p_opts', 'E_opts');
