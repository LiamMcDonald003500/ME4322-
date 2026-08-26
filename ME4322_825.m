% six bar linkage
% static equalibrium

clc;
clear;

% define joints
A = [7 4 0];
B = [5 16 0];
C = [25 25 0];
D = [23 10 0];
E = [18 35 0];
F = [43 32 0];
G = [45 17 0];

% LENGTH OF BARS
lAB = norm(B - A);
lBC = norm(C - B);
lCD = norm(D - C);
lBE = norm(E - B);
lEF = norm(F - E);
lFG = norm(G - F);

WAB = [0 -1 0];
WBEC = [0 -1 0];
WCD = [0 -1 0];
WEF = [0 -1 0];
WFG = [0 -1 0];

S1 = (A+B)/2;
S2 = (B+C+E)/3;
S3 = (C+D)/2;
S4 = (E+F)/2;
S5 = (F+G)/2;

syms FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin

ForceA = [FAx FAy 0];
ForceB = [FBx FBy 0];
ForceC = [FCx FCy 0];
ForceD = [FDx FDy 0];
ForceE = [FEx FEy 0];
ForceF = [FFx FFy 0];
ForceG = [FGx FGy 0];
InputTorque = [0 0 Tin];

AppliedForce = [50 0 0];

eqn1 = ForceA + ForceB + WAB == 0;

eqn2 = cross(A-S1,ForceA) + cross(B-S1,ForceB) + InputTorque == 0;

eqn3 = -ForceB + ForceC + ForceE + WBEC == 0;

eqn4 = cross(B-S2, -ForceB) + cross(C-S2, ForceC) + cross(E-S2, ForceE) == 0;

eqn5 = -ForceC + ForceD + WCD == 0;

eqn6 = cross(C-S3, -ForceC) + cross(D-S3, ForceD) == 0;

eqn7 = -ForceE + ForceF + WEF == 0;

eqn8 = cross(E-S4, -ForceE) + cross(F-S4, ForceF) == 0;

eqn9 = -ForceF + ForceG + WFG + AppliedForce == 0;

eqn10 = cross(F-S5, -ForceF) + cross(G-S5, ForceG) == 0;

%SOLVE THEM ALL!!!

eqnMatrix = [eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10];

StaticSolution=solve(eqnMatrix,[FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin]);

Force_Ay = double(StaticSolution.FAy);
Force_Ax = double(StaticSolution.FAx);
Force_By = double(StaticSolution.FBy);
Force_Bx = double(StaticSolution.FBx);
Force_Cy = double(StaticSolution.FCy);
Force_Cx = double(StaticSolution.FCx);
Force_Dy = double(StaticSolution.FDy);
Force_Dx = double(StaticSolution.FDx);
Force_Ey = double(StaticSolution.FEy);
Force_Ex = double(StaticSolution.FEx);
Force_Fy = double(StaticSolution.FFy);
Force_Fx = double(StaticSolution.FFx);
Force_Gy = double(StaticSolution.FGy);
Force_Gx = double(StaticSolution.FGx);

syms wBEC wCD wEF wFG

omegaAB = [0 0 1];
omegaBEC = [0 0 wBEC];
omegaCD = [0 0 wCD];


eqn11 = cross(omegaAB,B-A) + cross(omegaBEC,C-B) + cross(omegaCD,D-C) ==0;

loop1Solution = solve(eqn11,[wBEC wCD]);

angularVelocityBEC = double(loop1Solution.wBEC);
angularVelocityCD = double(loop1Solution.wCD);

omegaBEC2 = [0 0 angularVelocityBEC];
omegaCD2 = [0 0 angularVelocityCD];

omegaEF = [0 0 wEF];
omegaFG = [0 0 wFG];

eqn12 = cross(omegaCD2,C-D) + cross(omegaBEC2,E-C) + cross(omegaEF,F-E) + cross(omegaFG,G-F) == 0;

loop2Solution = solve(eqn12,[wEF wFG]);

angularVelocityEF = double(loop2Solution.wEF);
angularVelocityFG = double(loop2Solution.wFG);

syms aBEC aCD aFG aEF
alphaAB = [0 0 0];
alphaBEC = [0 0 aBEC];
alphaCD = [0 0 aCD];

a_BA = cross(alphaAB,B-A) + cross(omegaAB,cross(omegaAB,B-A));
a_CB = cross(alphaBEC,C-B) + cross(omegaBEC2,cross(omegaBEC2,C-B));
a_DC = cross(alphaCD,D-C) + cross(omegaCD2,cross(omegaCD2,D-C));

eqn13 = a_BA + a_CB + a_DC == 0;

loop1Accelu = solve(eqn13,[aBEC aCD]);

loop1AccelBEC = double(loop1Accelu.aBEC);
loop1AccelCD = double(loop1Accelu.aCD);

alphaBEC_vector = [0 0 loop1AccelBEC];
alphaCD_vector = [0 0 loop1AccelCD];

alphaEF = [0 0 aEF];
alphaFG = [0 0 aFG];

a_CD = cross(alphaCD_vector,C-D) + cross(omegaCD2,cross(omegaCD2,C-D));
a_EC = cross(alphaBEC_vector,E-C) + cross(omegaBEC2,cross(omegaBEC2,E-C));

angVelEF = [0 0 angularVelocityEF];
angVelFG = [0 0 angularVelocityFG];

a_FE = cross(alphaEF,F-E) + cross(angVelEF,cross(angVelEF,F-E));
a_GF = cross(alphaFG,G-F) + cross(angVelFG,cross(angVelFG,G-F));

eqn14 = a_CD + a_EC + a_FE + a_GF == 0;

loop2Accelu = solve(eqn14,[aEF aFG]);

%coding errs above please check and advise

vB_A = cross(omegaAB,B-A);

vE_B = cross(omegaBEC2,E-B);

vE_A = vE_B + vB_A;

%VELOCITY AT COM

V_S4_F = cross(angVelEF,S4-F);

V_F_G = cross(angVelFG,F-G);

vS4_G = V_S4_F + V_F_G;

%repeat for s1-5

