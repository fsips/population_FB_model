function v = fluxes_MIC(t,s,p,c, fixed_X)

% For numerical stability, don't continue with negative states
% s(s<0)        = 0; 
if nargin == 5
    s(find(~isnan(fixed_X))) = fixed_X(find(~isnan(fixed_X)));
end

%% States
v.V1          = s(1);       % Liter
v.V2          = s(2);
v.V3          = s(3);
v.V4          = s(4);
v.V5          = s(5);

v.B1          = s(6);       % 10^12 bacteria
v.B2          = s(7);
v.B3          = s(8);
v.B4          = s(9);
v.B5          = s(10);

v.F1          = s(11);      % 10^12 bacteria
v.F2          = s(12);
v.F3          = s(13);
v.F4          = s(14);
v.F5          = s(15);

v.SCFA_1      = s(16);      % mol
v.SCFA_2      = s(17);
v.SCFA_3      = s(18);
v.SCFA_4      = s(19);
v.SCFA_5      = s(20);

v.C_1        = s(21);
v.C_2        = s(22);
v.C_3        = s(23);
v.C_4        = s(24);
v.C_5        = s(25);

v.RS_1        = s(26);
v.RS_2        = s(27);
v.RS_3        = s(28);
v.RS_4        = s(29);
v.RS_5        = s(30);

v.PC_1        = s(31);
v.PC_2        = s(32);
v.PC_3        = s(33);
v.PC_4        = s(34);
v.PC_5        = s(35);

v.PU_1        = s(36);
v.PU_2        = s(37);
v.PU_3        = s(38);
v.PU_4        = s(39);
v.PU_5        = s(40);

v.SC_1        = s(41);
v.SC_2        = s(42);
v.SC_3        = s(43);
v.SC_4        = s(44);
v.SC_5        = s(45);

v.SU_1        = s(46);
v.SU_2        = s(47);
v.SU_3        = s(48);
v.SU_4        = s(49);
v.SU_5        = s(50);

v.d1          = s(51);
v.d2          = s(52);
v.d3          = s(53);
v.d4          = s(54);
v.d5          = s(55);

v.BA_pool     = s(56);

v.A1          = s(57);
v.A2          = s(58);
v.A3          = s(59);
v.A4          = s(60);
v.A5          = s(61);

v.cyc         = s(62);

%% Constants
v.Vin         = c.Vin;
v.SCFAin      = c.SCFAin;
v.Cin         = c.Cin;
v.RSin        = c.RSin;
v.Fin         = c.Fin;
v.Bin         = c.Bin;
v.k_co        = c.k_co;
v.k_uv        = c.k_uv;
v.k_ub        = c.k_ub;

v.tp          = c.tp;            
v.tc          = c.tc; 
v.rel_cI      = c.rel_cI;
v.rel_bO      = c.rel_bO;
v.pHopt       = c.pHopt;
v.pHinh       = c.pHinh;
v.k_inh_pH_B  = c.k_inh_pH_B;
v.k_inh_pH_F  = c.k_inh_pH_F;
v.SCFApC_B   = c.SCFApC_B;
v.SCFApC_F   = c.SCFApC_F;
v.Km          = c.Km;

v.raise_pH    = c.raise_pH;
v.Ka          = c.Ka;
v.BAopt       = c.BAopt;
v.BAinh       = c.BAinh;
v.k_inh_BA_B  = c.k_inh_BA_B;
v.FXR         = c.FXR;
v.cyc0        = c.cyc0 ;
v.dis         = c.dis;

v.ku          = max(c.ku-v.FXR*(v.cyc-v.cyc0),0) ;

%% Parameters

% Maximal growth
v.mu_maxB       = p(1);
v.mu_maxF       = p(2);

% Yield
v.Y_f           = p(3);
v.Y_b           = p(4);

% Carbohydrate fluxes
v.recyc_F       = p(5);
v.SCFAyield     = p(6);

% carbohydrate uptake
v.k_uS          = 0;
v.k_uA          = p(7);

% Bile acids
v.p_dec         = p(8);
v.p_deh         = p(9);

% Buffer and deacidification
v.b             = p(10);
v.frac_d        = p(11);

% Bile acid toxicity
v.inh_fact_F    = p(12);

% Firmicutes substrate preference
v.c_factor_F    = p(13);

%% Transluminal uptake of carbohydrates

% Uptake of SCFA
v.uSCFA_1   = v.k_uA * v.SCFA_1;
v.uSCFA_2   = v.k_uA * v.SCFA_2;
v.uSCFA_3   = v.k_uA * v.SCFA_3;
v.uSCFA_4   = v.k_uA * v.SCFA_4;
v.uSCFA_5   = v.k_uA * v.SCFA_5;

% Uptake of C
v.uC_1     = v.k_uS * v.C_1;
v.uC_2     = v.k_uS * v.C_2;
v.uC_3     = v.k_uS * v.C_3;
v.uC_4     = v.k_uS * v.C_4;
v.uC_5     = v.k_uS * v.C_5;

%% Bile acid kinetics

% Calculate BA inhibition of F
v.k_inh_BA_F    = v.k_inh_BA_B / v.inh_fact_F;              % BA inhibition of F

% Deconjugation
v.P1_dec       = v.p_dec *(v.B1+v.F1) * v.PC_1;
v.P2_dec       = v.p_dec *(v.B2+v.F2) * v.PC_2;
v.P3_dec       = v.p_dec *(v.B3+v.F3) * v.PC_3;
v.P4_dec       = v.p_dec *(v.B4+v.F4) * v.PC_4;
v.P5_dec       = v.p_dec *(v.B5+v.F5) * v.PC_5;

v.S1_dec       = v.p_dec *(v.B1+v.F1) * v.SC_1;
v.S2_dec       = v.p_dec *(v.B2+v.F2) * v.SC_2;
v.S3_dec       = v.p_dec *(v.B3+v.F3) * v.SC_3;
v.S4_dec       = v.p_dec *(v.B4+v.F4) * v.SC_4;
v.S5_dec       = v.p_dec *(v.B5+v.F5) * v.SC_5;

% Transformation
v.BA1_deh      = v.p_deh  *(v.F1)* v.PU_1;
v.BA2_deh      = v.p_deh  *(v.F2)* v.PU_2;
v.BA3_deh      = v.p_deh  *(v.F3)* v.PU_3;
v.BA4_deh      = v.p_deh  *(v.F4)* v.PU_4;
v.BA5_deh      = v.p_deh  *(v.F5)* v.PU_5;

% Passive uptake
v.P1_pu       = v.k_ub * v.PU_1;
v.P2_pu       = v.k_ub * v.PU_2;
v.P3_pu       = v.k_ub * v.PU_3;
v.P4_pu       = v.k_ub * v.PU_4;
v.P5_pu       = v.k_ub * v.PU_5;

v.S1_pu       = v.k_ub * v.SU_1;
v.S2_pu       = v.k_ub * v.SU_2;
v.S3_pu       = v.k_ub * v.SU_3;
v.S4_pu       = v.k_ub * v.SU_4;
v.S5_pu       = v.k_ub * v.SU_5;

% Input
v.cI          = v.rel_cI*v.rel_bO*v.BA_pool;
v.cO          = (v.S1_pu+v.S2_pu+v.S3_pu+v.S4_pu+v.S5_pu + v.P1_pu+v.P2_pu+v.P3_pu+v.P4_pu+v.P5_pu);
v.fsc         = (v.S1_pu+v.S2_pu+v.S3_pu+v.S4_pu+v.S5_pu)/v.cO;

term1         = (v.fsc * v.cO) / ((1-v.rel_cI)*(v.rel_bO*v.BA_pool + v.ku) + v.cO);
term2         = 1 - ((1-v.rel_cI)*v.rel_bO*v.BA_pool) /  ((1-v.rel_cI)*(v.rel_bO*v.BA_pool + v.ku ) + v.cO);
v.fS          = term1 / term2;

% v.fS        = (v.rel_cI*v.cO/v.cI*v.fsc)/(1-((v.rel_bO*v.BA_pool-v.ku)/(v.rel_bO*v.BA_pool)));

v.iP          = v.rel_cI*v.ku + v.rel_cI*v.rel_bO*v.BA_pool*(1-v.fS);
v.iS          = v.rel_cI*v.rel_bO*v.BA_pool*v.fS;


% Calculation of relative BA level (in mmol/L)
v.relBA_1      = ((v.tp*v.PU_1+v.tc*v.PC_1) + (v.SU_1+v.tc*v.SC_1))/v.V1/1e3;
v.relBA_2      = ((v.tp*v.PU_2+v.tc*v.PC_2) + (v.SU_2+v.tc*v.SC_2))/v.V2/1e3;
v.relBA_3      = ((v.tp*v.PU_3+v.tc*v.PC_3) + (v.SU_3+v.tc*v.SC_3))/v.V3/1e3;
v.relBA_4      = ((v.tp*v.PU_4+v.tc*v.PC_4) + (v.SU_4+v.tc*v.SC_4))/v.V4/1e3;
v.relBA_5      = ((v.tp*v.PU_5+v.tc*v.PC_5) + (v.SU_5+v.tc*v.SC_5))/v.V5/1e3;

%% pH kinetics
% Volume uptake
v.uV_1      = v.k_uv * v.V1;
v.uV_2      = v.k_uv * v.V2;
v.uV_3      = v.k_uv * v.V3;
v.uV_4      = v.k_uv * v.V4;
v.uV_5      = v.k_uv * v.V5;

% Input of buffer acid and d 
v.p_A1      = v.b;                              % Input of A- + HA, mol/min
v.x_in      = (v.p_A1) /(1+1/(10^(7+v.raise_pH-(-1*log10(v.Ka)))));   
v.id        = v.x_in - (10^-(7+v.raise_pH));

v.x_1      = (-1*(v.Ka-v.d1/v.V1)+sqrt(((v.Ka-v.d1/v.V1)^2)+4*(v.SCFA_1+v.A1)/v.V1*v.Ka))/2;
v.x_2      = (-1*(v.Ka-v.d2/v.V2)+sqrt(((v.Ka-v.d2/v.V2)^2)+4*(v.SCFA_2+v.A2)/v.V2*v.Ka))/2;
v.x_3      = (-1*(v.Ka-v.d3/v.V3)+sqrt(((v.Ka-v.d3/v.V3)^2)+4*(v.SCFA_3+v.A3)/v.V3*v.Ka))/2;
v.x_4      = (-1*(v.Ka-v.d4/v.V4)+sqrt(((v.Ka-v.d4/v.V4)^2)+4*(v.SCFA_4+v.A4)/v.V4*v.Ka))/2;
v.x_5      = (-1*(v.Ka-v.d5/v.V5)+sqrt(((v.Ka-v.d5/v.V5)^2)+4*(v.SCFA_5+v.A5)/v.V5*v.Ka))/2;

v.dd1       = v.uSCFA_1 * v.frac_d;
v.dd2       = v.uSCFA_2 * v.frac_d;
v.dd3       = v.uSCFA_3 * v.frac_d;
v.dd4       = v.uSCFA_4 * v.frac_d;
v.dd5       = v.uSCFA_5 * v.frac_d;

% Calculation of pH
v.pH1       = -1*log10(max(v.x_1 - v.d1/v.V1, 10^-7));
v.pH2       = -1*log10(max(v.x_2 - v.d2/v.V2, 10^-7));
v.pH3       = -1*log10(max(v.x_3 - v.d3/v.V3, 10^-7));
v.pH4       = -1*log10(max(v.x_4 - v.d4/v.V4, 10^-7));
v.pH5       = -1*log10(max(v.x_5 - v.d5/v.V5, 10^-7));

%% Microbiota kinetics

% Microbial growth: B
%   pH inhibition factors
v.in_B1_pH    = 1./(1+((v.pH1-v.pHopt))^v.pHinh*v.k_inh_pH_B);
v.in_B2_pH    = 1./(1+((v.pH2-v.pHopt))^v.pHinh*v.k_inh_pH_B);
v.in_B3_pH    = 1./(1+((v.pH3-v.pHopt))^v.pHinh*v.k_inh_pH_B);
v.in_B4_pH    = 1./(1+((v.pH4-v.pHopt))^v.pHinh*v.k_inh_pH_B);
v.in_B5_pH    = 1./(1+((v.pH5-v.pHopt))^v.pHinh*v.k_inh_pH_B);

%   BA inhibition factors
if log10(v.relBA_1)>v.BAopt,  v.in_B1_BA    = 1./(1+(log10(v.relBA_1)-v.BAopt).^v.BAinh*v.k_inh_BA_B); else v.in_B1_BA    = 1; end
if log10(v.relBA_2)>v.BAopt,  v.in_B2_BA    = 1./(1+(log10(v.relBA_2)-v.BAopt).^v.BAinh*v.k_inh_BA_B); else v.in_B2_BA    = 1; end
if log10(v.relBA_3)>v.BAopt,  v.in_B3_BA    = 1./(1+(log10(v.relBA_3)-v.BAopt).^v.BAinh*v.k_inh_BA_B); else v.in_B3_BA    = 1; end
if log10(v.relBA_4)>v.BAopt,  v.in_B4_BA    = 1./(1+(log10(v.relBA_4)-v.BAopt).^v.BAinh*v.k_inh_BA_B); else v.in_B4_BA    = 1; end
if log10(v.relBA_5)>v.BAopt,  v.in_B5_BA    = 1./(1+(log10(v.relBA_5)-v.BAopt).^v.BAinh*v.k_inh_BA_B); else v.in_B5_BA    = 1; end

%   Growth rates
v.mu_B1       = v.mu_maxB * (v.C_1/v.V1) / (v.Km + (v.C_1/v.V1)) * v.in_B1_pH * v.in_B1_BA;
v.mu_B2       = v.mu_maxB * (v.C_2/v.V2) / (v.Km + (v.C_2/v.V2)) * v.in_B2_pH * v.in_B2_BA;
v.mu_B3       = v.mu_maxB * (v.C_3/v.V3) / (v.Km + (v.C_3/v.V3)) * v.in_B3_pH * v.in_B3_BA;
v.mu_B4       = v.mu_maxB * (v.C_4/v.V4) / (v.Km + (v.C_4/v.V4)) * v.in_B4_pH * v.in_B4_BA;
v.mu_B5       = v.mu_maxB * (v.C_5/v.V5) / (v.Km + (v.C_5/v.V5)) * v.in_B5_pH * v.in_B5_BA;

%   Proliferation
v.p_B1        = v.mu_B1 * v.B1;
v.p_B2        = v.mu_B2 * v.B2;
v.p_B3        = 0; %v.mu_B3 * v.B3;
v.p_B4        = 0; %v.mu_B4 * v.B4;
v.p_B5        = 0; %v.mu_B5 * v.B5;

% Microbial growth: F
%   pH inhibition factors
v.in_F1_pH    = 1./(1+((v.pH1-v.pHopt))^v.pHinh*v.k_inh_pH_F);
v.in_F2_pH    = 1./(1+((v.pH2-v.pHopt))^v.pHinh*v.k_inh_pH_F);
v.in_F3_pH    = 1./(1+((v.pH3-v.pHopt))^v.pHinh*v.k_inh_pH_F);
v.in_F4_pH    = 1./(1+((v.pH4-v.pHopt))^v.pHinh*v.k_inh_pH_F);
v.in_F5_pH    = 1./(1+((v.pH5-v.pHopt))^v.pHinh*v.k_inh_pH_F);

%   BA inhibition factors
if log10(v.relBA_1)>v.BAopt,  v.in_F1_BA    = 1./(1+(log10(v.relBA_1)-v.BAopt).^v.BAinh*v.k_inh_BA_F); else v.in_F1_BA    = 1; end
if log10(v.relBA_2)>v.BAopt,  v.in_F2_BA    = 1./(1+(log10(v.relBA_2)-v.BAopt).^v.BAinh*v.k_inh_BA_F); else v.in_F2_BA    = 1; end
if log10(v.relBA_3)>v.BAopt,  v.in_F3_BA    = 1./(1+(log10(v.relBA_3)-v.BAopt).^v.BAinh*v.k_inh_BA_F); else v.in_F3_BA    = 1; end
if log10(v.relBA_4)>v.BAopt,  v.in_F4_BA    = 1./(1+(log10(v.relBA_4)-v.BAopt).^v.BAinh*v.k_inh_BA_F); else v.in_F4_BA    = 1; end
if log10(v.relBA_5)>v.BAopt,  v.in_F5_BA    = 1./(1+(log10(v.relBA_5)-v.BAopt).^v.BAinh*v.k_inh_BA_F); else v.in_F5_BA    = 1; end

%   Growth rates
v.mu_F1       = v.mu_maxF * ((v.c_factor_F*v.C_1+v.RS_1)/v.V1) / (v.Km + ((v.c_factor_F*v.C_1+v.RS_1)/v.V1)) * v.in_F1_pH * v.in_F1_BA;
v.mu_F2       = v.mu_maxF * ((v.c_factor_F*v.C_2+v.RS_2)/v.V2) / (v.Km + ((v.c_factor_F*v.C_2+v.RS_2)/v.V2)) * v.in_F2_pH * v.in_F2_BA;
v.mu_F3       = v.mu_maxF * ((v.c_factor_F*v.C_3+v.RS_3)/v.V3) / (v.Km + ((v.c_factor_F*v.C_3+v.RS_3)/v.V3)) * v.in_F3_pH * v.in_F3_BA;
v.mu_F4       = v.mu_maxF * ((v.c_factor_F*v.C_4+v.RS_4)/v.V4) / (v.Km + ((v.c_factor_F*v.C_4+v.RS_4)/v.V4)) * v.in_F4_pH * v.in_F4_BA;
v.mu_F5       = v.mu_maxF * ((v.c_factor_F*v.C_5+v.RS_5)/v.V5) / (v.Km + ((v.c_factor_F*v.C_5+v.RS_5)/v.V5)) * v.in_F5_pH * v.in_F5_BA;

%   Proliferation
v.p_F1        = v.mu_F1 * v.F1;
v.p_F2        = v.mu_F2 * v.F2;
v.p_F3        = 0; %v.mu_F3 * v.F3;
v.p_F4        = 0; %v.mu_F4 * v.F4;
v.p_F5        = 0; %v.mu_F5 * v.F5;

%% Carbohydrate metabolism

% RS depletion
%   Currently, B does not break down RS
if (v.C_1+v.RS_1)>0, v.dRS1_dt_F = (v.mu_F1 * (v.RS_1)./(v.c_factor_F*v.C_1+v.RS_1)) * v.F1 / v.Y_f; else v.dRS1_dt_F = 0; end
if (v.C_2+v.RS_2)>0, v.dRS2_dt_F = (v.mu_F2 * (v.RS_2)./(v.c_factor_F*v.C_2+v.RS_2)) * v.F2 / v.Y_f; else v.dRS2_dt_F = 0; end
if (v.C_3+v.RS_3)>0, v.dRS3_dt_F = (v.mu_F3 * (v.RS_3)./(v.c_factor_F*v.C_3+v.RS_3)) * v.F3 / v.Y_f; else v.dRS3_dt_F = 0; end
if (v.C_4+v.RS_4)>0, v.dRS4_dt_F = (v.mu_F4 * (v.RS_4)./(v.c_factor_F*v.C_4+v.RS_4)) * v.F4 / v.Y_f; else v.dRS4_dt_F = 0; end
if (v.C_5+v.RS_5)>0, v.dRS5_dt_F = (v.mu_F5 * (v.RS_5)./(v.c_factor_F*v.C_5+v.RS_5)) * v.F5 / v.Y_f; else v.dRS5_dt_F = 0; end

% C Fermentation
if (v.C_1+v.RS_1)>0, v.dC1_dt_F = (v.mu_F1 * (v.c_factor_F*v.C_1)./(v.c_factor_F*v.C_1+v.RS_1)) * v.F1 / v.Y_f;   else v.dC1_dt_F = 0;  end
if (v.C_2+v.RS_2)>0, v.dC2_dt_F = (v.mu_F2 * (v.c_factor_F*v.C_2)./(v.c_factor_F*v.C_2+v.RS_2)) * v.F2 / v.Y_f;   else v.dC2_dt_F = 0;  end
if (v.C_3+v.RS_3)>0, v.dC3_dt_F = (v.mu_F3 * (v.c_factor_F*v.C_3)./(v.c_factor_F*v.C_3+v.RS_3)) * v.F3 / v.Y_f;   else v.dC3_dt_F = 0;  end
if (v.C_4+v.RS_4)>0, v.dC4_dt_F = (v.mu_F4 * (v.c_factor_F*v.C_4)./(v.c_factor_F*v.C_4+v.RS_4)) * v.F4 / v.Y_f;   else v.dC4_dt_F = 0;  end
if (v.C_5+v.RS_5)>0, v.dC5_dt_F = (v.mu_F5 * (v.c_factor_F*v.C_5)./(v.c_factor_F*v.C_5+v.RS_5)) * v.F5 / v.Y_f;   else v.dC5_dt_F = 0;  end

v.dC1_dt_B = v.mu_B1 * v.B1 / v.Y_b;
v.dC2_dt_B = v.mu_B2 * v.B2 / v.Y_b;
v.dC3_dt_B = v.mu_B3 * v.B3 / v.Y_b;
v.dC4_dt_B = v.mu_B4 * v.B4 / v.Y_b;
v.dC5_dt_B = v.mu_B5 * v.B5 / v.Y_b;


%% Transit
v.tr_B_2_1    = v.k_co * v.B1;
v.tr_B_3_2    = v.k_co * v.B2;
v.tr_B_4_3    = v.k_co * v.B3;
v.tr_B_5_4    = v.k_co * v.B4;
v.tr_B_6_5    = v.k_co * v.B5;

v.tr_F_2_1    = v.k_co * v.F1;
v.tr_F_3_2    = v.k_co * v.F2;
v.tr_F_4_3    = v.k_co * v.F3;
v.tr_F_5_4    = v.k_co * v.F4;
v.tr_F_6_5    = v.k_co * v.F5;

v.tr_SCFA_2_1 = v.k_co * v.SCFA_1;
v.tr_SCFA_3_2 = v.k_co * v.SCFA_2;
v.tr_SCFA_4_3 = v.k_co * v.SCFA_3;
v.tr_SCFA_5_4 = v.k_co * v.SCFA_4;
v.tr_SCFA_6_5 = v.k_co * v.SCFA_5;

v.tr_C_2_1   = v.k_co * v.C_1;
v.tr_C_3_2   = v.k_co * v.C_2;
v.tr_C_4_3   = v.k_co * v.C_3;
v.tr_C_5_4   = v.k_co * v.C_4;
v.tr_C_6_5   = v.k_co * v.C_5;

v.tr_RS_2_1   = v.k_co * v.RS_1;
v.tr_RS_3_2   = v.k_co * v.RS_2;
v.tr_RS_4_3   = v.k_co * v.RS_3;
v.tr_RS_5_4   = v.k_co * v.RS_4;
v.tr_RS_6_5   = v.k_co * v.RS_5;

v.tr_V_2_1    = v.k_co * v.V1;
v.tr_V_3_2    = v.k_co * v.V2;
v.tr_V_4_3    = v.k_co * v.V3;
v.tr_V_5_4    = v.k_co * v.V4;
v.tr_V_6_5    = v.k_co * v.V5;

v.tr_d_2_1    = v.k_co * v.d1;
v.tr_d_3_2    = v.k_co * v.d2;
v.tr_d_4_3    = v.k_co * v.d3;
v.tr_d_5_4    = v.k_co * v.d4;
v.tr_d_6_5    = v.k_co * v.d5;

v.tr_PU_2_1    = v.k_co * v.PU_1;
v.tr_PU_3_2    = v.k_co * v.PU_2;
v.tr_PU_4_3    = v.k_co * v.PU_3;
v.tr_PU_5_4    = v.k_co * v.PU_4;
v.tr_PU_6_5    = v.k_co * v.PU_5;

v.tr_PC_2_1    = v.k_co * v.PC_1;
v.tr_PC_3_2    = v.k_co * v.PC_2;
v.tr_PC_4_3    = v.k_co * v.PC_3;
v.tr_PC_5_4    = v.k_co * v.PC_4;
v.tr_PC_6_5    = v.k_co * v.PC_5;

v.tr_SU_2_1    = v.k_co * v.SU_1;
v.tr_SU_3_2    = v.k_co * v.SU_2;
v.tr_SU_4_3    = v.k_co * v.SU_3;
v.tr_SU_5_4    = v.k_co * v.SU_4;
v.tr_SU_6_5    = v.k_co * v.SU_5;

v.tr_SC_2_1    = v.k_co * v.SC_1;
v.tr_SC_3_2    = v.k_co * v.SC_2;
v.tr_SC_4_3    = v.k_co * v.SC_3;
v.tr_SC_5_4    = v.k_co * v.SC_4;
v.tr_SC_6_5    = v.k_co * v.SC_5;

v.tr_A_2_1    = v.k_co * v.A1;
v.tr_A_3_2    = v.k_co * v.A2;
v.tr_A_4_3    = v.k_co * v.A3;
v.tr_A_5_4    = v.k_co * v.A4;
v.tr_A_6_5    = v.k_co * v.A5;

%% SUMS

v.synthesis    = max(v.ku-v.FXR*(v.cyc-v.cyc0),0);
v.synthesis_dif= -v.FXR*(v.cyc-v.cyc0);

v.SCFA         = v.SCFA_1 + v.SCFA_2 + v.SCFA_3 + v.SCFA_4 + v.SCFA_5;
v.C            = v.C_1 + v.C_2 + v.C_3 + v.C_4 + v.C_5;
v.RS           = v.RS_1 + v.RS_2 + v.RS_3 + v.RS_4 + v.RS_5;

v.B            = v.B1 + v.B2 + v.B3 + v.B4 + v.B5;
v.F            = v.F1 + v.F2 + v.F3 + v.F4 + v.F5;

v.PU           = v.PU_1 + v.PU_2 + v.PU_3 + v.PU_4 + v.PU_5;
v.PC           = v.PC_1 + v.PC_2 + v.PC_3 + v.PC_4 + v.PC_5;
v.SU           = v.SU_1 + v.SU_2 + v.SU_3 + v.SU_4 + v.SU_5;
v.SC           = v.SC_1 + v.SC_2 + v.SC_3 + v.SC_4 + v.SC_5;

v.P            = v.PU + v.PC;
v.S            = v.SU + v.SC;
v.cOP          = (v.P1_pu+v.P2_pu+v.P3_pu+v.P4_pu+v.P5_pu);
v.cOS          = (v.S1_pu+v.S2_pu+v.S3_pu+v.S4_pu+v.S5_pu);
