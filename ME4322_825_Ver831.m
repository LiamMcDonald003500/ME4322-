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

syms aBEC aCD aFG aEF



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


alphaAB = [0 0 0];
alphaBEC = [0 0 aBEC];
alphaCD = [0 0 aCD];
alphaEF = [0 0 aEF];
alphaFG = [0 0 aFG];

a_BA = cross(alphaAB,B-A) + cross(omegaAB,cross(omegaAB,B-A));
a_CB = cross(alphaBEC,C-B) + cross(omegaBEC2,cross(omegaBEC2,C-B));
a_DC = cross(alphaCD,D-C) + cross(omegaCD2,cross(omegaCD2,D-C));

eqn13 = a_BA + a_CB + a_DC == 0;

loop1Accelu = solve(eqn13,[aBEC aCD]);

loop1AccelBEC = double(loop1Accelu.aBEC);
loop1AccelCD = double(loop1Accelu.aCD);

alphaBEC_vector = [0 0 loop1AccelBEC];
alphaCD_vector = [0 0 loop1AccelCD];


a_CD = cross(alphaCD_vector,C-D) + cross(omegaCD2,cross(omegaCD2,C-D));
a_EC = cross(alphaBEC_vector,E-C) + cross(omegaBEC2,cross(omegaBEC2,E-C));

angVelEF = [0 0 angularVelocityEF];
angVelFG = [0 0 angularVelocityFG];

a_FE = cross(alphaEF,F-E) + cross(angVelEF,cross(angVelEF,F-E));
a_FG = cross(alphaFG,G-F) + cross(angVelFG,cross(angVelFG,G-F));

eqn14 = a_CD + a_EC + a_FE + a_FG == 0;

loop2Accelu = solve(eqn14,[aEF aFG]);

AEF = double(loop2Accelu.aEF);
AFG = double(loop2Accelu.aFG);


%temp corr last two

vB_A = cross(omegaAB,B-A);

vE_B = cross(omegaBEC2,E-B);

vE_A = vE_B + vB_A;

%VELOCITY AT COM

V_S4_F = cross(angVelEF,S4-F);

V_F_G = cross(angVelFG,F-G);

vS4_G = V_S4_F + V_F_G;

aS1_A = cross(alphaAB,S1-A) + cross(omegaAB,cross(omegaAB,S1-A));

aS2_A = cross(alphaBEC_vector,S2-B) + cross(omegaBEC2,cross(omegaBEC2,S2-B)) + a_BA;

aS3_D = cross(alphaCD_vector,S3-D) + cross(omegaCD2,cross(omegaCD2,S3-D));

aS4_F = cross([0 0 AEF],S4-F) + cross(angVelEF,cross(angVelEF,S4-F)) +cross([0 0 AFG],F-G) + cross(angVelFG,cross(angVelFG,F-G));

aS5_G = cross([0 0 AFG],S5-G) + cross(angVelFG,cross(angVelFG,S5-G));

%MMOI

J_AB = 1;
J_BEC = 1;
J_CD = 1;
J_EF = 1;
J_FG = 1;

MassAB = 1;
MassBEC = 1;
MassCD = 1;
MassEF = 1;
MassFG = 1;

%overhead may be incomplete
syms NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy NFEx NFEy NFFx NFFy NFGx NFGy NTin
NForceA = [NFAx NFAy 0];
NForceB = [NFBx NFBy 0];
NForceC = [NFCx NFCy 0];
NForceD = [NFDx NFDy 0];
NForceE = [NFEx NFEy 0];
NForceF = [NFFx NFFy 0];
NForceG = [NFGx NFGy 0];
NinputTorque = [0 0 NTin];

%link AB
eqn15 = NForceA + NForceB + WAB == MassAB * aS1_A;

%sigmaM=0
eqn16 = cross(A-S1,NForceA) + cross(B-S1, NForceB) + NinputTorque == J_AB * alphaAB;

%link BEC
eqn17 = -NForceB + NForceC + NForceE + WBEC == MassBEC * aS2_A;

eqn18 = cross(B-S2,-NForceB) + cross(C-S2,NForceC) + cross(E-S2,NForceE) == J_BEC * alphaBEC_vector;

eqn19 = -NForceC + NForceD + WCD == MassCD * aS3_D;

eqn20 = cross(C-S3,-NForceC) + cross(D-S3,NForceD) == J_CD * alphaCD_vector;

eqn21 = -NForceE + NForceF + WEF == MassEF * aS4_F;

eqn22 = cross(E-S4,-NForceE) + cross(F-S4,NForceF) == J_EF * [0 0 AEF];

eqn23 = -NForceF + NForceG + WFG + AppliedForce == MassFG * aS5_G;

eqn24 = cross(F-S5,-NForceF) + cross(G-S5,NForceG) == J_FG * [0 0 AFG];

solnmatrix = [eqn15,eqn16,eqn17,eqn18,eqn19,eqn20,eqn21,eqn22,eqn23,eqn24];
dynamicsolution = solve(solnmatrix, [NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy NFEx NFEy NFFx NFFy NFGx NFGy NTin]);

NForceAx = double(dynamicsolution.NFAx);
NForceAy = double(dynamicsolution.NFAy);
NForceBx = double(dynamicsolution.NFBx);
NForceBy = double(dynamicsolution.NFBy);
NForceCx = double(dynamicsolution.NFCx);
NForceCy = double(dynamicsolution.NFCy);
NForceDx = double(dynamicsolution.NFDx);
NForceDy = double(dynamicsolution.NFDy);
NForceEx = double(dynamicsolution.NFEx);
NForceEy = double(dynamicsolution.NFEy);
NForceFx = double(dynamicsolution.NFFx);
NForceFy = double(dynamicsolution.NFFy);
NForceGx = double(dynamicsolution.NFGx);
NForceGy = double(dynamicsolution.NFGy);
NTorque = double(dynamicsolution.NTin);

StaticForceA = [Force_Ax,Force_Ay];
StaticForceB = [Force_Bx,Force_By];
StaticForceC = [Force_Cx,Force_Cy];
StaticForceD = [Force_Dx,Force_Dy];
StaticForceE = [Force_Ex,Force_Ey];
StaticForceF = [Force_Fx,Force_Fy];
StaticForceG = [Force_Gx,Force_Gy];


DynamicForceA = [NForceAx, NForceAy];
DynamicForceB = [NForceBx, NForceBy];
DynamicForceC = [NForceCx, NForceCy];
DynamicForceD = [NForceDx, NForceDy];
DynamicForceE = [NForceEx, NForceEy];
DynamicForceF = [NForceFx, NForceFy];
DynamicForceG = [NForceGx, NForceGy];

StaticMag = [
    norm(StaticForceA)
    norm(StaticForceB)
    norm(StaticForceC)
    norm(StaticForceD)
    norm(StaticForceE)
    norm(StaticForceF)
    norm(StaticForceG)];

DynamicMag = [
    norm(DynamicForceA)
    norm(DynamicForceB)
    norm(DynamicForceC)
    norm(DynamicForceD)
    norm(DynamicForceE)
    norm(DynamicForceF)
    norm(DynamicForceG)];

ForceDiff = DynamicMag - StaticMag;

Joint = ["A"; "B"; "C"; "D"; "E"; "F"; "G"];

ForceComp = table(Joint,StaticMag,DynamicMag,ForceDiff,'VariableNames',{'Joint','Static_Force','Dynamic_Force','Difference'});

disp(' ')
disp('FORCE COMPARISON ---------------------------------')
disp(ForceComp)

%Circle intersection

%joint coords already defined as [joint name]
%lengths defined as l[link name]

initial_theta = atan2(B(2)-A(2),B(1)-A(1));

if (initial_theta<0)
    inputanglecir = 2*pi + initial_theta;
else
    inputanglecir = initial_theta;
end

for theta=1:1:360
    %New joint B

    b_new = [lAB*cos(inputanglecir+deg2rad(theta)) lAB*sin(inputanglecir+deg2rad(theta)) 0];
    
    %new position of C, b_new as ctr and BC as radius, then D as center and
    %DC as radius

    [Cx, Cy] = circcirc(b_new(1),b_new(2),lBC,D(1),D(2),lCD);
    
end

%concept steps: convert points to polar, find initial angle, use towards a
%point_new statement before using circle analysis for (new, new, length,
%point, length)
%prof plz check and advise (maybe extention for midnight reqd)

