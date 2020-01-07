function [E, s, d] = output_MIC(v)

%% Microbiota outputs
% 1. total microbiota in feces: 2-5*10^11 per gr wet weight
%    -> average of 2.5*10^11 per gr wet weight
%    -> of which 90 % is B+F
%    -> wet weight: 100 grams / day
%    -> is 100/(60*24) = 0.0694 grams / min
%    -> so per minute we expect
%    2.5e11 * (100/(60*24)) * 0.9 B+F
d(1)    = (2.5e11 * (100/(60*24)) * 0.9)/1e12;
s(1)    = (v.tr_B_6_5 + v.tr_F_6_5) ./ (v.tr_V_6_5*1e3);

% 2. F/B in feces is 1
d(2)    = 1;
if v.tr_B_6_5 > 0
    s(2)    = v.tr_F_6_5/v.tr_B_6_5;
    if s(2) > 1e6
        s(2)= 1e6;
    end
else
    s(2)    = 1e6;
end

%% pH outputs
% 3. pH in the cecum
d(3)    = 10.^(-1*5.6);
s(3)    = 10.^(-1*v.pH1);

% 4. Terminal pH
d(4)    = 10.^(-1*6.6);
s(4)    = 10.^(-1*v.pH5);

%% Bile acids
% 5. Primary BA in feces
s(5)    = (v.tr_PU_6_5 + v.tr_PC_6_5)/(v.tr_PU_6_5 + v.tr_PC_6_5 + v.tr_SU_6_5 + v.tr_SC_6_5);
d(5)    = 0.05;

% 6. Conjugated BA in feces
s(6)    = (v.tr_SC_6_5 + v.tr_PC_6_5)/(v.tr_PU_6_5 + v.tr_PC_6_5 + v.tr_SU_6_5 + v.tr_SC_6_5);
d(6)    = 0.02;

% 7. Recycled BA
s(7)    = (v.tr_PU_6_5 + v.tr_PC_6_5 + v.tr_SU_6_5 + v.tr_SC_6_5) / (v.iP + v.iS);
d(7)    = 1-0.28;

%% SCFA
% SCFA concentration in co1
s(8)    = v.SCFA_1 / v.V1;
d(8)    = 100/1e3;

% SCFA concentration in co5
s(9)    = v.SCFA_5 / v.V5;
d(9)    = 40/1e3;

% SCFA in feces
s(10)   = v.tr_SCFA_6_5;
d(10)   = 10/1e3   / (60*24);


%%
E       = [(s-d)./d];
E       = [E(:)];


