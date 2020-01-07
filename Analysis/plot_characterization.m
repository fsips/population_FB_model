function [v_SS] = plot_characterization(p,c,x0)
%PLOT_CHARACTERIZATION Plot steady state model results and residuals
%   Inputs:
%   p       - parameter vector
%   c       - constants 
%   x0      - initial conditions
%
%   Outputs:
%   v_SS    - variables in steady state

%% Run simulation
global bt
bt = tic;   % For time limit on ode - simulation will not run without

options = odeset('NonNegative',1, 'AbsTol', 1e-9, 'RelTol', 1e-6);
[t,s] = ode15s(@ode_MIC,[0 1e3*24*60],x0,options,p,c,[]);

%% Plot 1: Show states to verify model has reached steady state
figure();
subplot(2,1,1)
plot(t,s);

for it = size(s,1)
    v(it) = fluxes_MIC(t(it), s(it,:),p,c);
end

%% Plot 2: Show inhibition by pH and BA
h = figure();
ys = 1;
xs = 2;

c_B     = [1/5 4/5 1];
c_F     = [2/5 4/5 7/10];
c_C1    = [2/5 4/5  3/10];
c_C2    = [3/5 9/10 5/10];
c_C3    = [4/5 1    7/10];
cm_C    = [c_C1; c_C2; c_C3];

subplot(ys, xs, 1)
pH  = 4:0.1:7;
inh_B = 1./(1+((pH-7)./1).^4*v(end).k_inh_pH_B);
inh_F = 1./(1+((pH-7)./1).^4*v(end).k_inh_pH_F);
plot(pH, inh_B, 'Color', c_B, 'LineWidth', 2); hold on
plot(pH, inh_F, 'Color', c_F, 'LineWidth', 2);
legend({'B' 'F'}, 'Location', 'SouthEast')
xlim([4 7])
xlabel('pH')
ylabel('Relative growth rate')
axis square

subplot(ys, xs, 2)
rBA  = 10.^(-3:0.1:1);
inh_B = 1./(1+(log10(rBA)-v(end).BAopt).^4*v(end).k_inh_BA_B);
inh_F = 1./(1+(log10(rBA)-v(end).BAopt).^4*v(end).k_inh_BA_F);
semilogx(rBA, inh_B, 'Color', c_B, 'LineWidth', 2); hold on
semilogx(rBA, inh_F, 'Color', c_F, 'LineWidth', 2);
xlim(10.^[-3 1])
xlabel('[BA] (mmol/L)')
ylabel('Relative growth rate')
set(gca, 'XTick', 10.^[-3 -2 -1 0 1])
set(gca, 'XTickLabel', 10.^[-3 -2 -1 0 1])
axis square

%% Plot 3: Characetrize what is happening in the colon
h = figure();
ys = 2;
xs = 3;

% A. Total bacterial content
subplot(ys, xs, 1)
bar([1 2 3 4 5], [v(end).F1 v(end).F2 v(end).F3 v(end).F4 v(end).F5]+[v(end).B1 v(end).B2 v(end).B3 v(end).B4 v(end).B5], 'FaceColor', [59 191 199]/255)
ylabel('Microbiota (# . 10^{12})')
set(gca, 'XTick', [1 2 3 4 5])
set(gca, 'XTickLabel', {'co1' 'co2' 'co3' 'co4' 'co5'})
axis square

% B. SCFA content
subplot(ys, xs, 2)
bar([1 2 3 4 5], [v(end).SCFA_1 v(end).SCFA_2 v(end).SCFA_3 v(end).SCFA_4 v(end).SCFA_5]')
legend({'SCFA' 'C' 'RS'})
ylabel('Carbohydrates (molar eq)')
set(gca, 'XTick', [1 2 3 4 5])
ylim([0 0.03])
colormap(cm_C)
set(gca, 'XTickLabel', {'co1' 'co2' 'co3' 'co4' 'co5'})
axis square

% C. pH
subplot(ys, xs, 3)
bar([1 2 3 4 5], [v(end).pH1 v(end).pH2 v(end).pH3 v(end).pH4 v(end).pH5], 'FaceColor', [85 88 165]/255)
ylabel('pH')
set(gca, 'XTick', [1 2 3 4 5])
set(gca, 'XTickLabel', {'co1' 'co2' 'co3' 'co4' 'co5'})
axis square

% D. F/B
subplot(ys, xs, 4)
bar([1 2 3 4 5 6], [v(end).F1 v(end).F2 v(end).F3 v(end).F4 v(end).F5 v(end).tr_F_6_5]./[v(end).B1 v(end).B2 v(end).B3 v(end).B4 v(end).B5 v(end).tr_B_6_5], 'FaceColor', [59 191 199]/255)
ylabel('F/B')
set(gca, 'XTick', [1 2 3 4 5 6])
set(gca, 'XTickLabel', {'co1' 'co2' 'co3' 'co4' 'co5' 'fe'})
axis square

% E. Bile acid content
subplot(ys, xs, 5)
BA = [ v(end).PU_1 v(end).PU_2  v(end).PU_3  v(end).PU_4  v(end).PU_5
    v(end).PC_1 v(end).PC_2  v(end).PC_3  v(end).PC_4  v(end).PC_5
    v(end).SU_1 v(end).SU_2  v(end).SU_3  v(end).SU_4  v(end).SU_5
    v(end).SC_1 v(end).SC_2  v(end).SC_3  v(end).SC_4  v(end).SC_5]';
h = area([1 2 3 4 5], BA);
set(h(1),'FaceColor',[190 186 205]/255)
set(h(2),'FaceColor',[146 137 204]/255)
set(h(3),'FaceColor',[190 172 193]/255)
set(h(4),'FaceColor',[175 113 188]/255)

legend({'PU' 'PC' 'SU' 'SC'})
ylabel('BA (\mumol)')
set(gca, 'XTick', [1 2 3 4 5])
set(gca, 'XTickLabel', {'co1' 'co2' 'co3' 'co4' 'co5'})
axis square


%% F. Residuals
[E, s, d] = output_MIC(v(end));
figure()
plot([1:10],(d-s)./d, 'sk', 'MarkerFaceColor', [0 0 0]); hold on
xlabel('Datapoint')
ylabel('(d_i - s_i) / d_i')
sum(E.^2)
v_SS = v(end);