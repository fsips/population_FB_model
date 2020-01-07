function dX = ode_MIC(t, s, p, c, fixed_X)

if nargin<5
    fixed_X = ones(size(s))*NaN;
end

% For optimization, cut simulation short if it takes too long
global bt       
et = toc(bt);
if et > 240*10
    error('Simulation is taking too long')
    t
end

% For numerical stability, don't continue with negative states
% s(s<0)        = 0; 

% Calculate right side
v = fluxes_MIC(t,s,p,c, fixed_X);

%% Differential equations: V
dV1_dt      = v.Vin  - v.uV_1              - v.tr_V_2_1;
dV2_dt      =        - v.uV_2 + v.tr_V_2_1 - v.tr_V_3_2;
dV3_dt      =        - v.uV_3 + v.tr_V_3_2 - v.tr_V_4_3;
dV4_dt      =        - v.uV_4 + v.tr_V_4_3 - v.tr_V_5_4;
dV5_dt      =        - v.uV_5 + v.tr_V_5_4 - v.tr_V_6_5;

%% Differential equations: H
dd1_dt      = v.id      + v.dd1               - v.tr_d_2_1;
dd2_dt      =           + v.dd2  + v.tr_d_2_1 - v.tr_d_3_2;
dd3_dt      =           + v.dd3  + v.tr_d_3_2 - v.tr_d_4_3;
dd4_dt      =           + v.dd4  + v.tr_d_4_3 - v.tr_d_5_4;
dd5_dt      =           + v.dd5  + v.tr_d_5_4 - v.tr_d_6_5;

dA1dt       = + v.p_A1 - v.tr_A_2_1             ;
dA2dt       =          - v.tr_A_3_2 + v.tr_A_2_1;
dA3dt       =          - v.tr_A_4_3 + v.tr_A_3_2;
dA4dt       =          - v.tr_A_5_4 + v.tr_A_4_3;
dA5dt       =          - v.tr_A_6_5 + v.tr_A_5_4;

%% Differential equations: Bacteriodetes
dB1dt       = + v.p_B1 - v.tr_B_2_1 + v.Bin     ;
dB2dt       = + v.p_B2 - v.tr_B_3_2 + v.tr_B_2_1;
dB3dt       = + v.p_B3 - v.tr_B_4_3 + v.tr_B_3_2;
dB4dt       = + v.p_B4 - v.tr_B_5_4 + v.tr_B_4_3;
dB5dt       = + v.p_B5 - v.tr_B_6_5 + v.tr_B_5_4;

%% Differential equations: Firmicutes
dF1dt       = + v.p_F1 - v.tr_F_2_1 + v.Fin     ;
dF2dt       = + v.p_F2 - v.tr_F_3_2 + v.tr_F_2_1;
dF3dt       = + v.p_F3 - v.tr_F_4_3 + v.tr_F_3_2;
dF4dt       = + v.p_F4 - v.tr_F_5_4 + v.tr_F_4_3;
dF5dt       = + v.p_F5 - v.tr_F_6_5 + v.tr_F_5_4;

%% Differential equations: Carbohydrates
dSCFA_1_dt  = v.SCFAin              - v.tr_SCFA_2_1 + v.SCFApC_B*v.dC1_dt_B*v.SCFAyield +  v.SCFApC_F*(v.dC1_dt_F+v.dRS1_dt_F)*v.SCFAyield - v.uSCFA_1;
dSCFA_2_dt  =       + v.tr_SCFA_2_1 - v.tr_SCFA_3_2 + v.SCFApC_B*v.dC2_dt_B*v.SCFAyield +  v.SCFApC_F*(v.dC2_dt_F+v.dRS2_dt_F)*v.SCFAyield - v.uSCFA_2;
dSCFA_3_dt  =       + v.tr_SCFA_3_2 - v.tr_SCFA_4_3 + v.SCFApC_B*v.dC3_dt_B*v.SCFAyield +  v.SCFApC_F*(v.dC3_dt_F+v.dRS3_dt_F)*v.SCFAyield - v.uSCFA_3;
dSCFA_4_dt  =       + v.tr_SCFA_4_3 - v.tr_SCFA_5_4 + v.SCFApC_B*v.dC4_dt_B*v.SCFAyield +  v.SCFApC_F*(v.dC4_dt_F+v.dRS4_dt_F)*v.SCFAyield - v.uSCFA_4;
dSCFA_5_dt  =       + v.tr_SCFA_5_4 - v.tr_SCFA_6_5 + v.SCFApC_B*v.dC5_dt_B*v.SCFAyield +  v.SCFApC_F*(v.dC5_dt_F+v.dRS5_dt_F)*v.SCFAyield - v.uSCFA_5;

dC_1_dt    = v.Cin              - v.tr_C_2_1 - v.dC1_dt_B - v.dC1_dt_F + v.dRS1_dt_F*v.recyc_F - v.uC_1;
dC_2_dt    =       + v.tr_C_2_1 - v.tr_C_3_2 - v.dC2_dt_B - v.dC2_dt_F + v.dRS2_dt_F*v.recyc_F - v.uC_2;
dC_3_dt    =       + v.tr_C_3_2 - v.tr_C_4_3 - v.dC3_dt_B - v.dC3_dt_F + v.dRS3_dt_F*v.recyc_F - v.uC_3;
dC_4_dt    =       + v.tr_C_4_3 - v.tr_C_5_4 - v.dC4_dt_B - v.dC4_dt_F + v.dRS4_dt_F*v.recyc_F - v.uC_4;
dC_5_dt    =       + v.tr_C_5_4 - v.tr_C_6_5 - v.dC5_dt_B - v.dC5_dt_F + v.dRS5_dt_F*v.recyc_F - v.uC_5;

dRS_1_dt    = v.RSin              - v.tr_RS_2_1 - v.dRS1_dt_F*(1+v.recyc_F);
dRS_2_dt    =       + v.tr_RS_2_1 - v.tr_RS_3_2 - v.dRS2_dt_F*(1+v.recyc_F);
dRS_3_dt    =       + v.tr_RS_3_2 - v.tr_RS_4_3 - v.dRS3_dt_F*(1+v.recyc_F);
dRS_4_dt    =       + v.tr_RS_4_3 - v.tr_RS_5_4 - v.dRS4_dt_F*(1+v.recyc_F);
dRS_5_dt    =       + v.tr_RS_5_4 - v.tr_RS_6_5 - v.dRS5_dt_F*(1+v.recyc_F);

%% Differential equations: Bile acids
dPC_1_dt    = + v.iP                - v.tr_PC_2_1                 - v.P1_dec;
dPC_2_dt    =                       - v.tr_PC_3_2 + v.tr_PC_2_1   - v.P2_dec;
dPC_3_dt    =                       - v.tr_PC_4_3 + v.tr_PC_3_2   - v.P3_dec;
dPC_4_dt    =                       - v.tr_PC_5_4 + v.tr_PC_4_3   - v.P4_dec;
dPC_5_dt    =                       - v.tr_PC_6_5 + v.tr_PC_5_4   - v.P5_dec;

dPU_1_dt    =                       - v.tr_PU_2_1                 + v.P1_dec - v.BA1_deh - v.P1_pu;
dPU_2_dt    =                       - v.tr_PU_3_2 + v.tr_PU_2_1   + v.P2_dec - v.BA2_deh - v.P2_pu;
dPU_3_dt    =                       - v.tr_PU_4_3 + v.tr_PU_3_2   + v.P3_dec - v.BA3_deh - v.P3_pu;
dPU_4_dt    =                       - v.tr_PU_5_4 + v.tr_PU_4_3   + v.P4_dec - v.BA4_deh - v.P4_pu;
dPU_5_dt    =                       - v.tr_PU_6_5 + v.tr_PU_5_4   + v.P5_dec - v.BA5_deh - v.P5_pu;

dSC_1_dt    = + v.iS                - v.tr_SC_2_1                 - v.S1_dec;
dSC_2_dt    =                       - v.tr_SC_3_2 + v.tr_SC_2_1   - v.S2_dec;
dSC_3_dt    =                       - v.tr_SC_4_3 + v.tr_SC_3_2   - v.S3_dec;
dSC_4_dt    =                       - v.tr_SC_5_4 + v.tr_SC_4_3   - v.S4_dec;
dSC_5_dt    =                       - v.tr_SC_6_5 + v.tr_SC_5_4   - v.S5_dec;

dSU_1_dt    =                       - v.tr_SU_2_1                 + v.S1_dec + v.BA1_deh - v.S1_pu;
dSU_2_dt    =                       - v.tr_SU_3_2 + v.tr_SU_2_1   + v.S2_dec + v.BA2_deh - v.S2_pu;
dSU_3_dt    =                       - v.tr_SU_4_3 + v.tr_SU_3_2   + v.S3_dec + v.BA3_deh - v.S3_pu;
dSU_4_dt    =                       - v.tr_SU_5_4 + v.tr_SU_4_3   + v.S4_dec + v.BA4_deh - v.S4_pu;
dSU_5_dt    =                       - v.tr_SU_6_5 + v.tr_SU_5_4   + v.S5_dec + v.BA5_deh - v.S5_pu;

d_BApool_dt = v.ku                  - v.tr_PC_6_5 - v.tr_PU_6_5 - v.tr_SC_6_5 - v.tr_SU_6_5 - v.dis*v.BA_pool;

d_cyc       =   1/100 * (((1-v.rel_cI)*v.rel_bO*v.BA_pool*((1-v.fS)*10+v.fS)+(1-v.rel_cI)*v.ku*10) - v.cyc);

%% Combine to vector
dX          = [dV1_dt;dV2_dt;dV3_dt;dV4_dt;dV5_dt;dB1dt;dB2dt;dB3dt;dB4dt;dB5dt;                % 1-10
    dF1dt;dF2dt;dF3dt;dF4dt;dF5dt;dSCFA_1_dt;dSCFA_2_dt;dSCFA_3_dt;dSCFA_4_dt;dSCFA_5_dt;       % 11-20
    dC_1_dt;dC_2_dt;dC_3_dt;dC_4_dt;dC_5_dt;dRS_1_dt;dRS_2_dt;dRS_3_dt;dRS_4_dt;dRS_5_dt        % 21-30
    dPC_1_dt;dPC_2_dt;dPC_3_dt;dPC_4_dt;dPC_5_dt;dPU_1_dt;dPU_2_dt;dPU_3_dt;dPU_4_dt;dPU_5_dt;  % 31-40
    dSC_1_dt;dSC_2_dt;dSC_3_dt;dSC_4_dt;dSC_5_dt;dSU_1_dt;dSU_2_dt;dSU_3_dt;dSU_4_dt;dSU_5_dt;  % 41-50
    dd1_dt;dd2_dt;dd3_dt;dd4_dt;dd5_dt;d_BApool_dt;dA1dt;dA2dt;dA3dt;dA4dt;dA5dt;d_cyc];              % 51-61


