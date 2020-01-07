function [p] = parameters_MIC()

% Microbiota kinetics: maximal growth rates
mu_maxB         = 1/60; % Min 0
mu_maxF         = 1/120;% Min 0

% Microbiota kinetics: yield coefficients
Y_f             = 1;    % Min 0
Y_b             = 50;    % Min 0

% Carbohydrate metabolism: RS recycling; SCFA yield
recyc_F         = 0.1;    % Min 0
SCFAyield       = 0.2;    % Min 0, Max 1
p_SCFA_uptake   = 10^(-2.6252);

% Bile acid metabolism
p_deconjugation = 10^(-3.4255)*5;
p_dehydroxylat  = 10^(-3.1431)*5;

% pH modifier
p_pH            = 1e2/1e3/(24*60);   % Buffer in (mol/day)
p_pH2           = 0.1       ;        % Fraction of SCFA uptake for which HCO3- is exchanged

% BA inhibition factor
inh_fact_F      = (1/16)/ 0.01;

% Substrate preference factor for C
c_factor_F      = 0.01;

p = [mu_maxB
     mu_maxF
     Y_f 
     Y_b
     recyc_F
     SCFAyield
     p_SCFA_uptake
     p_deconjugation  % Deconjugation
     p_dehydroxylat   % Dehydroxylation (by F)
     p_pH
     p_pH2   
     inh_fact_F
     c_factor_F
     ];
