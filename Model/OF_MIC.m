function E = OF_MIC(p0, c, x0)

p1 = p0;

global bt       
bt = tic;

try
    options = odeset('NonNegative',1, 'AbsTol', 1e-9, 'RelTol', 1e-6);
    
    % Initial simulation, with fixed carbohydrate states
    fixed_X = fixed_states(x0);
    [t1,s1] = ode15s(@ode_MIC,[0 1e3*24*60],x0,options,p1,c,fixed_X);
    
    % Second simulation, with initial conditions distilled from first 
    x02     =  s1(end,:);
    x02(find(~isnan(fixed_X))) = fixed_X(find(~isnan(fixed_X)));
    [t2,s2] = ode15s(@ode_MIC,[0 1e3*24*60],x02,options,p1,c);
    
    if t2(end) ~= 1e3*24*60
        error()
    end
    
    for it = size(s2,1)
        v2(it) = fluxes_MIC(t2(it), s2(it,:),p1,c);
    end
    
    E       = output_MIC(v2(end));
    E_pad   = zeros(20,1);
    E       = [E(:); E_pad];
catch
    E       = ones(30,1)*1e5;
    E       = E(:);
end