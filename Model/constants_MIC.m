function [c] = constants_MIC()

%% Inputs

c.Vin       = 1      / (60*24); % Volume input           (L/min)
c.SCFAin    = 0      / (60*24); % SCFA input             (mol / min)
c.Cin       = 0.075  / (60*24); % Carbohydrate input     (mol / min)
c.RSin      = 0.025  / (60*24); % Resistant starch input (mol / min)
c.Fin       = 0.9*0.5*2.5e9/1e12/(60*24); % Bacterial input (#10^12 bacteria / min)
c.Bin       = 0.9*0.5*2.5e9/1e12/(60*24); % Bacterial input (#10^12 bacteria / min)

%% Transits and uptake

c.k_co      = 2.12e-3;          % Colon transit time     (1/min)
c.k_uv      = 1.24e-3;          % Volume uptake rate     (1/min)
c.k_ub      = 1.92e-4;          % Bile acid uptake       (1/min)

%% Bile acid metabolism

c.ku        = 0.82;             % Bile acid synthesis    (umol / min)
c.rel_cI    = 0.045;            % To colon of secreted
c.rel_bO    = 6 / (60*24);      % Biliary output (rel to pool size)

%% Inhibition of microbial proliferation

c.tc        = 1/100;            % Decreased toxicity of conjugated BA vs uncojugated
c.tp        = 1/10;             % Decreased toxicity of (unconjugated) primary BA vs secondary

c.BAopt     = -3;               % Top of BA inhibition curve
c.BAinh     = 4;                % Power of BA inhibition curve
c.k_inh_BA_B= 1/16;             % BA inhibition of B

c.pHopt     = 7;                % Top of PH inhibition curve              
c.pHinh     = 4;                % Power of PH inhibition curve 
c.k_inh_pH_B= 0.625;            % pH inhibition of B
c.k_inh_pH_F= 0.125;            % pH inhibition of F

%% Carbohydrate metabolism

c.SCFApC_B  = 1.25;             % SCFA multiplication factor (rel to F) (B)              
c.SCFApC_F  = 1;                % SCFA multiplication factor (F) 
c.Km        = 50e-6;            % Km (50 muM = 50e-6 M)

%% Buffer constants

c.raise_pH  = -0.5;
c.Ka        = 10^(-4.8);

%% FXR regulation

c.FXR       = 0;
c.cyc0      = 0;
c.dis       = 0;
