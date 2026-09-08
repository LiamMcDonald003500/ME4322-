clc;
clear;
close all;

%% Joint coordinates
A0 = [1.4 0.485];
B0 = [1.67 0.99];
C0 = [0.255 1.035];
D0 = [0.285 0.055];
E0 = [0.195 2.54];
F0 = [-0.98 2.57];
G0 = [0.05 0.2];

% Grounded joints
A = A0.';
D = D0.';
G = G0.';

%% Link lengths
lAB = norm(B0-A0);
lBC = norm(C0-B0);
lDC = norm(C0-D0);
lDE = norm(E0-D0);
lEF = norm(F0-E0);
lFG = norm(G0-F0);

% Gripper location: 1.843 m beyond F along the E-to-F direction
gripperDistance = 1.843;
gripper0 = F0 + gripperDistance*(F0-E0)/lEF;

%% Link masses
% Each link is 10 cm wide and 5 mm thick, and each joint hole is 6 cm
% in diameter. The assumed material is AISI 1095 steel.

function mass = linkmass(L)

% Link dimensions are in meters.
width = 0.1;
thickness = 0.005;
jointdiameter = 0.06;
density = 7850000; % g/m^3

V_rec = width*thickness*L;
V_hole = 2*pi*((jointdiameter/2)^2)*thickness;
V_net = V_rec - V_hole;

mass = (density*V_net)/1000;

end

% All masses are in kilograms.
mAB = linkmass(lAB);
mBC = linkmass(lBC);
mDC = linkmass(lDC);
mDE = linkmass(lDE);
mEF = linkmass(lEF);
mFG = linkmass(lFG);

%% Mass moments of inertia
I_AB = (1/12)*mAB*lAB^2;
I_BC = (1/12)*mBC*lBC^2;
I_DC = (1/12)*mDC*lDC^2;
I_DE = (1/12)*mDE*lDE^2;
I_EF = (1/12)*mEF*lEF^2;
I_FG = (1/12)*mFG*lFG^2;

% Combined properties of rigid link CDE
mCDE = mDC + mDE;
S_DC = (D0 + C0)/2;
S_DE = (D0 + E0)/2;
S_CDE = (mDC*S_DC + mDE*S_DE)/mCDE;
I_CDE = I_DC + mDC*norm(S_DC-S_CDE)^2 + I_DE + mDE*norm(S_DE-S_CDE)^2;

%% Initial angles
thetaAB0 = atan2(B0(2)-A0(2), B0(1)-A0(1));
thetaDC = atan2(C0(2)-D0(2), C0(1)-D0(1));
thetaDE = atan2(E0(2)-D0(2), E0(1)-D0(1));

% Constant angle between DC and DE
betaCDE = thetaDE - thetaDC;

%% Input rotation
inputchange = deg2rad(0:1:360);
thetaAB = thetaAB0 + inputchange;

numbersteps = length(thetaAB);

% Input speed required for 12,500 cycles in 9 hours
omega_required = (12500*2*pi)/(9*60*60);
omegaAB = omega_required; % rad/s
time = inputchange/omegaAB;

% Joint-coordinate arrays
B = NaN(2, numbersteps);
C = NaN(2, numbersteps);
E = NaN(2, numbersteps);
F = NaN(2, numbersteps);

validPosition = false(1, numbersteps);

previousC = C0.';
previousF = F0.';

for k = 1:numbersteps

    B(:,k) = A + lAB*[cos(thetaAB(k)); sin(thetaAB(k))];

    [xC,yC] = circcirc(B(1,k),B(2,k),lBC,D(1),D(2),lDC);

    C1 = [xC(1); yC(1)];
    C2 = [xC(2); yC(2)];

    if any(isnan(xC)) || any(isnan(yC))
        fprintf("No valid position at input change %.1f degrees.\n",rad2deg(inputchange(k)));
        continue
    end

    distanceC1 = norm(C1 - previousC);
    distanceC2 = norm(C2 - previousC);

    if distanceC1 <= distanceC2
        C(:,k) = C1;
    else
        C(:,k) = C2;
    end

    thetaDC = atan2(C(2,k)-D(2), C(1,k)-D(1));
    thetaDE = thetaDC + betaCDE;

    E(:,k) = D + lDE*[cos(thetaDE); sin(thetaDE)];

    [xF,yF] = circcirc(E(1,k),E(2,k),lEF,G(1),G(2),lFG);

    F1 = [xF(1); yF(1)];
    F2 = [xF(2); yF(2)];

    if any(isnan(xF)) || any(isnan(yF))
        fprintf("Point F cannot be assembled at %.1f degrees.\n",rad2deg(inputchange(k)));
        B(:,k) = NaN;
        C(:,k) = NaN;
        E(:,k) = NaN;
        continue
    end

    distanceF1 = norm(F1 - previousF);
    distanceF2 = norm(F2 - previousF);

    if distanceF1 <= distanceF2
        F(:,k) = F1;
    else
        F(:,k) = F2;
    end

    previousC = C(:,k);
    previousF = F(:,k);

    validPosition(k) = true;
end

if ~all(validPosition)
    warning("Some input positions are invalid. Nearby velocities and accelerations will be NaN.");
end

%% Joint velocity analysis

dt = time(2) - time(1);

vB = [gradient(B(1,:),dt);
      gradient(B(2,:),dt)];

vC = [gradient(C(1,:),dt);
      gradient(C(2,:),dt)];

vE = [gradient(E(1,:),dt);
      gradient(E(2,:),dt)];

vF = [gradient(F(1,:),dt);
      gradient(F(2,:),dt)];

% Grounded joints
vA = zeros(2,numbersteps);
vD = zeros(2,numbersteps);
vG = zeros(2,numbersteps);

%% Joint acceleration analysis

aB = [gradient(vB(1,:),dt);
      gradient(vB(2,:),dt)];

aC = [gradient(vC(1,:),dt);
      gradient(vC(2,:),dt)];

aE = [gradient(vE(1,:),dt);
      gradient(vE(2,:),dt)];

aF = [gradient(vF(1,:),dt);
      gradient(vF(2,:),dt)];

% Grounded joints
aA = zeros(2,numbersteps);
aD = zeros(2,numbersteps);
aG = zeros(2,numbersteps);

% Magnitudes at every input position

speedA = vecnorm(vA,2,1);
speedB = vecnorm(vB,2,1);
speedC = vecnorm(vC,2,1);
speedD = vecnorm(vD,2,1);
speedE = vecnorm(vE,2,1);
speedF = vecnorm(vF,2,1);
speedG = vecnorm(vG,2,1);

accelA = vecnorm(aA,2,1);
accelB = vecnorm(aB,2,1);
accelC = vecnorm(aC,2,1);
accelD = vecnorm(aD,2,1);
accelE = vecnorm(aE,2,1);
accelF = vecnorm(aF,2,1);
accelG = vecnorm(aG,2,1);

%% Joint position plots

inputDegrees = rad2deg(inputchange);

figure;

subplot(2,1,1)
plot(inputDegrees,[B(1,:);C(1,:);E(1,:);F(1,:)].','LineWidth',1.4);
grid on;
xlabel('Input Rotation (degrees)');
ylabel('x Position (m)');
title('Joint x Positions');
legend('B','C','E','F','Location','best');

subplot(2,1,2)
plot(inputDegrees,[B(2,:);C(2,:);E(2,:);F(2,:)].','LineWidth',1.4);
grid on;
xlabel('Input Rotation (degrees)');
ylabel('y Position (m)');
title('Joint y Positions');
legend('B','C','E','F','Location','best');

saveas(gcf,'joint_positions.png');

%% Joint trajectory plot

figure;
plot(B(1,:),B(2,:),'LineWidth',1.4);
hold on;
plot(C(1,:),C(2,:),'LineWidth',1.4);
plot(E(1,:),E(2,:),'LineWidth',1.4);
plot(F(1,:),F(2,:),'LineWidth',1.4);

plot(A(1),A(2),'ko','MarkerFaceColor','k');
plot(D(1),D(2),'ko','MarkerFaceColor','k');
plot(G(1),G(2),'ko','MarkerFaceColor','k');

axis equal;
grid on;
xlabel('x Position (m)');
ylabel('y Position (m)');
title('Joint Trajectories');
legend('B','C','E','F','A','D','G','Location','bestoutside');

saveas(gcf,'joint_trajectories.png');

%% Joint velocity plot

figure;
jointSpeeds = [speedA;speedB;speedC;speedD;speedE;speedF;speedG];
plot(inputDegrees,jointSpeeds.','LineWidth',1.4);

grid on;
xlabel('Input Rotation (degrees)');
ylabel('Joint Speed (m/s)');
title('Joint Velocity Magnitudes');
legend('A','B','C','D','E','F','G','Location','bestoutside');

saveas(gcf,'joint_velocities.png');

%% Joint acceleration plot

figure;
jointAccelerations = [accelA;accelB;accelC;accelD;accelE;accelF;accelG];
plot(inputDegrees,jointAccelerations.','LineWidth',1.4);

grid on;
xlabel('Input Rotation (degrees)');
ylabel('Joint Acceleration (m/s^2)');
title('Joint Linear Acceleration Magnitudes');
legend('A','B','C','D','E','F','G','Location','bestoutside');

saveas(gcf,'joint_accelerations.png');

%% Link angular velocity and acceleration

% Link orientation angles at every position
angleAB = unwrap(atan2(B(2,:)-A(2),B(1,:)-A(1)));
angleBC = unwrap(atan2(C(2,:)-B(2,:),C(1,:)-B(1,:)));

% CDE is one rigid link, so line DC defines its orientation
angleCDE = unwrap(atan2(C(2,:)-D(2),C(1,:)-D(1)));
angleEF = unwrap(atan2(F(2,:)-E(2,:),F(1,:)-E(1,:)));
angleFG = unwrap(atan2(F(2,:)-G(2),F(1,:)-G(1)));

%% Link angular velocities

omegaAB_plot = gradient(angleAB,dt);
omegaBC = gradient(angleBC,dt);
omegaCDE = gradient(angleCDE,dt);
omegaEF = gradient(angleEF,dt);
omegaFG = gradient(angleFG,dt);

%% Link angular accelerations

alphaAB = gradient(omegaAB_plot,dt);
alphaBC = gradient(omegaBC,dt);
alphaCDE = gradient(omegaCDE,dt);
alphaEF = gradient(omegaEF,dt);
alphaFG = gradient(omegaFG,dt);

%% Angular velocity plot

figure;

linkAngularVelocities = [omegaAB_plot;omegaBC;omegaCDE;omegaEF;omegaFG];
plot(inputDegrees,linkAngularVelocities.','LineWidth',1.4);

grid on;
xlabel('Input Rotation (degrees)');
ylabel('Angular Velocity (rad/s)');
title('Link Angular Velocities');
legend('AB','BC','CDE','EF','FG','Location','bestoutside');

saveas(gcf,'link_angular_velocities.png');

%% Angular acceleration plot

figure;

linkAngularAccelerations = [alphaAB;alphaBC;alphaCDE;alphaEF;alphaFG];
plot(inputDegrees,linkAngularAccelerations.','LineWidth',1.4);

grid on;
xlabel('Input Rotation (degrees)');
ylabel('Angular Acceleration (rad/s^2)');
title('Link Angular Accelerations');
legend('AB','BC','CDE','EF','FG','Location','bestoutside');

saveas(gcf,'link_angular_accelerations.png');

%% Link center-of-mass positions

A_all = repmat(A,1,numbersteps);
D_all = repmat(D,1,numbersteps);
G_all = repmat(G,1,numbersteps);

S_AB_all = (A_all+B)/2;
S_BC_all = (B+C)/2;

% CDE is composed of segments DC and DE
S_DC_all = (D_all+C)/2;
S_DE_all = (D_all+E)/2;

S_CDE_all = (mDC*S_DC_all + mDE*S_DE_all)/mCDE;

S_EF_all = (E+F)/2;
S_FG_all = (F+G_all)/2;

% Gripper position throughout the cycle
gripper = F + gripperDistance*(F-E)/lEF;


%% Link center-of-mass velocities

vS_AB = [gradient(S_AB_all(1,:),dt);
         gradient(S_AB_all(2,:),dt)];

vS_BC = [gradient(S_BC_all(1,:),dt);
         gradient(S_BC_all(2,:),dt)];

vS_CDE = [gradient(S_CDE_all(1,:),dt);
          gradient(S_CDE_all(2,:),dt)];

vS_EF = [gradient(S_EF_all(1,:),dt);
         gradient(S_EF_all(2,:),dt)];

vS_FG = [gradient(S_FG_all(1,:),dt);
         gradient(S_FG_all(2,:),dt)];


%% Link center-of-mass accelerations

aS_AB = [gradient(vS_AB(1,:),dt);
         gradient(vS_AB(2,:),dt)];

aS_BC = [gradient(vS_BC(1,:),dt);
         gradient(vS_BC(2,:),dt)];

aS_CDE = [gradient(vS_CDE(1,:),dt);
          gradient(vS_CDE(2,:),dt)];

aS_EF = [gradient(vS_EF(1,:),dt);
         gradient(vS_EF(2,:),dt)];

aS_FG = [gradient(vS_FG(1,:),dt);
         gradient(vS_FG(2,:),dt)];

% Center-of-mass acceleration magnitudes
accelS_AB = vecnorm(aS_AB,2,1);
accelS_BC = vecnorm(aS_BC,2,1);
accelS_CDE = vecnorm(aS_CDE,2,1);
accelS_EF = vecnorm(aS_EF,2,1);
accelS_FG = vecnorm(aS_FG,2,1);

%% First-position static analysis

% Initial joint coordinates in 3D
A3 = [A0 0];
B3 = [B0 0];
C3 = [C0 0];
D3 = [D0 0];
E3 = [E0 0];
F3 = [F0 0];
G3 = [G0 0];
gripper3 = [gripper0 0];

% Link centers of mass in 3D
S_AB3 = [(A0+B0)/2 0];
S_BC3 = [(B0+C0)/2 0];
S_CDE3 = [S_CDE 0];
S_EF3 = [(E0+F0)/2 0];
S_FG3 = [(F0+G0)/2 0];

syms FAx FAy FBx FBy FCx FCy
syms FDx FDy FEx FEy FFx FFy FGx FGy Tin

ForceA = [FAx FAy 0];
ForceB = [FBx FBy 0];
ForceC = [FCx FCy 0];
ForceD = [FDx FDy 0];
ForceE = [FEx FEy 0];
ForceF = [FFx FFy 0];
ForceG = [FGx FGy 0];

InputTorque = [0 0 Tin];

% Weight vectors
g = 9.81;
GAB  = [0 -mAB*g 0];
GBC  = [0 -mBC*g 0];
GCDE = [0 -mCDE*g 0];
GEF  = [0 -mEF*g 0];
GFG  = [0 -mFG*g 0];

LoadGripper = [0 -200 0];

% Force equilibrium
staticForceAB  = ForceA + ForceB + GAB;
staticForceBC  = -ForceB + ForceC + GBC;
staticForceCDE = -ForceC + ForceD + ForceE + GCDE;
staticForceEF  = -ForceE + ForceF + GEF + LoadGripper;
staticForceFG  = -ForceF + ForceG + GFG;

% Moment equilibrium
% AB: moments about A
staticMomentAB = cross(B3-A3,ForceB) + cross(S_AB3-A3,GAB) + InputTorque;

% BC: moments about B
staticMomentBC = cross(C3-B3,ForceC) + cross(S_BC3-B3,GBC);

% CDE: moments about D
staticMomentCDE = cross(C3-D3,-ForceC) + cross(E3-D3,ForceE) + cross(S_CDE3-D3,GCDE);

% EF: moments about E
staticMomentEF = cross(F3-E3,ForceF) + cross(S_EF3-E3,GEF) + cross(gripper3-E3,LoadGripper);

% FG: moments about G
staticMomentFG = cross(F3-G3,-ForceF) + cross(S_FG3-G3,GFG);

staticEquations = [
    staticForceAB(1)  == 0
    staticForceAB(2)  == 0
    staticMomentAB(3) == 0

    staticForceBC(1)  == 0
    staticForceBC(2)  == 0
    staticMomentBC(3) == 0

    staticForceCDE(1)  == 0
    staticForceCDE(2)  == 0
    staticMomentCDE(3) == 0

    staticForceEF(1)  == 0
    staticForceEF(2)  == 0
    staticMomentEF(3) == 0

    staticForceFG(1)  == 0
    staticForceFG(2)  == 0
    staticMomentFG(3) == 0
];

unknowns = [FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin];

StaticSolution = solve(staticEquations,unknowns);

StaticForceValues = [
    double(StaticSolution.FAx)
    double(StaticSolution.FAy)
    double(StaticSolution.FBx)
    double(StaticSolution.FBy)
    double(StaticSolution.FCx)
    double(StaticSolution.FCy)
    double(StaticSolution.FDx)
    double(StaticSolution.FDy)
    double(StaticSolution.FEx)
    double(StaticSolution.FEy)
    double(StaticSolution.FFx)
    double(StaticSolution.FFy)
    double(StaticSolution.FGx)
    double(StaticSolution.FGy)
];

StaticInputTorque = double(StaticSolution.Tin);
%% First-position Newton's Second Law analysis

% Center-of-mass accelerations at the first position
aS_AB3   = [aS_AB(:,1).' 0];
aS_BC3   = [aS_BC(:,1).' 0];
aS_CDE3  = [aS_CDE(:,1).' 0];
aS_EF3   = [aS_EF(:,1).' 0];
aS_FG3   = [aS_FG(:,1).' 0];

% Angular accelerations at the first position
alphaAB3  = [0 0 alphaAB(1)];
alphaBC3  = [0 0 alphaBC(1)];
alphaCDE3 = [0 0 alphaCDE(1)];
alphaEF3  = [0 0 alphaEF(1)];
alphaFG3  = [0 0 alphaFG(1)];

% Newton's Second Law force equations
dynamicForceAB = ForceA + ForceB + GAB - mAB*aS_AB3;
dynamicForceBC = -ForceB + ForceC + GBC - mBC*aS_BC3;
dynamicForceCDE = -ForceC + ForceD + ForceE + GCDE - mCDE*aS_CDE3;
dynamicForceEF = -ForceE + ForceF + GEF + LoadGripper - mEF*aS_EF3;
dynamicForceFG = -ForceF + ForceG + GFG - mFG*aS_FG3;

% Newton-Euler moment equations about each mass center
dynamicMomentAB = cross(A3-S_AB3,ForceA) + cross(B3-S_AB3,ForceB) + InputTorque - I_AB*alphaAB3;
dynamicMomentBC = cross(B3-S_BC3,-ForceB) + cross(C3-S_BC3,ForceC) - I_BC*alphaBC3;
dynamicMomentCDE = cross(C3-S_CDE3,-ForceC) + cross(D3-S_CDE3,ForceD) + cross(E3-S_CDE3,ForceE) - I_CDE*alphaCDE3;
dynamicMomentEF = cross(E3-S_EF3,-ForceE) + cross(F3-S_EF3,ForceF) + cross(gripper3-S_EF3,LoadGripper) - I_EF*alphaEF3;
dynamicMomentFG = cross(F3-S_FG3,-ForceF) + cross(G3-S_FG3,ForceG) - I_FG*alphaFG3;

dynamicEquations = [
    dynamicForceAB(1)  == 0
    dynamicForceAB(2)  == 0
    dynamicMomentAB(3) == 0

    dynamicForceBC(1)  == 0
    dynamicForceBC(2)  == 0
    dynamicMomentBC(3) == 0

    dynamicForceCDE(1)  == 0
    dynamicForceCDE(2)  == 0
    dynamicMomentCDE(3) == 0

    dynamicForceEF(1)  == 0
    dynamicForceEF(2)  == 0
    dynamicMomentEF(3) == 0

    dynamicForceFG(1)  == 0
    dynamicForceFG(2)  == 0
    dynamicMomentFG(3) == 0
];

DynamicSolution = solve(dynamicEquations,unknowns);

DynamicForceValues = [
    double(DynamicSolution.FAx)
    double(DynamicSolution.FAy)
    double(DynamicSolution.FBx)
    double(DynamicSolution.FBy)
    double(DynamicSolution.FCx)
    double(DynamicSolution.FCy)
    double(DynamicSolution.FDx)
    double(DynamicSolution.FDy)
    double(DynamicSolution.FEx)
    double(DynamicSolution.FEy)
    double(DynamicSolution.FFx)
    double(DynamicSolution.FFy)
    double(DynamicSolution.FGx)
    double(DynamicSolution.FGy)
];

DynamicInputTorque = double(DynamicSolution.Tin);

%% First-position force and torque table

ResultName = [
    "F_Ax"
    "F_Ay"
    "F_Bx"
    "F_By"
    "F_Cx"
    "F_Cy"
    "F_Dx"
    "F_Dy"
    "F_Ex"
    "F_Ey"
    "F_Fx"
    "F_Fy"
    "F_Gx"
    "F_Gy"
    "Input torque"
];

Units = [
    repmat("N",14,1)
    "N*m"
];

StaticResults = [
    StaticForceValues
    StaticInputTorque
];

DynamicResults = [
    DynamicForceValues
    DynamicInputTorque
];

FirstPositionForceTable = table(ResultName,StaticResults,DynamicResults,Units)

%% Free-body diagrams

figure('Name','Link Free-Body Diagrams');
tiledlayout(2,3);

arrowLength = 0.18;

%% Link AB
nexttile;
plot([A0(1),B0(1)],[A0(2),B0(2)],'k-','LineWidth',5);
hold on;

drawForcePair(A0, 1,'A',arrowLength);
drawForcePair(B0, 1,'B',arrowLength);

S_AB = (A0+B0)/2;
drawWeight(S_AB,'W_{AB}',arrowLength);

% Symbolic input-torque label
text(S_AB(1),S_AB(2)+0.18,'\itT_{in} \circlearrowleft','FontSize',12,'HorizontalAlignment','center');

finishFBD('Link AB');


%% Link BC
nexttile;
plot([B0(1),C0(1)],[B0(2),C0(2)],'k-','LineWidth',5);
hold on;

% Force at B is opposite to the force shown on link AB
drawForcePair(B0,-1,'B',arrowLength);
drawForcePair(C0, 1,'C',arrowLength);

S_BC = (B0+C0)/2;
drawWeight(S_BC,'W_{BC}',arrowLength);

finishFBD('Link BC');


%% Rigid link CDE
nexttile;

% Outline of the rigid ternary link
plot([C0(1),D0(1),E0(1),C0(1)],[C0(2),D0(2),E0(2),C0(2)],'k-','LineWidth',5);
hold on;

drawForcePair(C0,-1,'C',arrowLength);
drawForcePair(D0, 1,'D',arrowLength);
drawForcePair(E0, 1,'E',arrowLength);

drawWeight(S_CDE,'W_{CDE}',arrowLength);

finishFBD('Rigid Link CDE');


%% Link EF
nexttile;
plot([E0(1),F0(1),gripper0(1)],[E0(2),F0(2),gripper0(2)],'k-','LineWidth',5);
hold on;

drawForcePair(E0,-1,'E',arrowLength);
drawForcePair(F0, 1,'F',arrowLength);

S_EF = (E0+F0)/2;
drawWeight(S_EF,'W_{EF}',arrowLength);

% External 200 N load applied at the gripper
quiver(gripper0(1),gripper0(2),0,-arrowLength,0,'Color',[0.85 0 0],'LineWidth',2,'MaxHeadSize',0.5);
text(gripper0(1)+0.03,gripper0(2)-arrowLength,'200 N','Color',[0.85 0 0],'FontWeight','bold');

finishFBD('Link EF');


%% Link FG
nexttile;
plot([F0(1),G0(1)],[F0(2),G0(2)],'k-','LineWidth',5);
hold on;

drawForcePair(F0,-1,'F',arrowLength);
drawForcePair(G0, 1,'G',arrowLength);

S_FG = (F0+G0)/2;
drawWeight(S_FG,'W_{FG}',arrowLength);

finishFBD('Link FG');


%% Remove unused sixth tile
nexttile;
axis off;

sgtitle('Free-Body Diagrams of Linkage');

saveas(gcf,'link_free_body_diagrams.png');

function drawForcePair(point,forceSign,jointName,L)

% Horizontal force component
quiver(point(1),point(2),forceSign*L,0,0,'b','LineWidth',1.5,'MaxHeadSize',0.5);

% Vertical force component
quiver(point(1),point(2),0,forceSign*L,0,'b','LineWidth',1.5,'MaxHeadSize',0.5);

if forceSign > 0
    xLabel = ['F_{',jointName,'x}'];
    yLabel = ['F_{',jointName,'y}'];
else
    xLabel = ['-F_{',jointName,'x}'];
    yLabel = ['-F_{',jointName,'y}'];
end

text(point(1)+forceSign*L,point(2)-0.04,xLabel,'Color','b','HorizontalAlignment','center');

text(point(1)+0.03,point(2)+forceSign*L,yLabel,'Color','b','VerticalAlignment','middle');

end


function drawWeight(center,weightName,L)

quiver(center(1),center(2),0,-L,0,'Color',[0.85 0 0],'LineWidth',1.5,'MaxHeadSize',0.5);

text(center(1)+0.03,center(2)-L,weightName,'Color',[0.85 0 0],'VerticalAlignment','top');

end


function finishFBD(linkTitle)

axis equal;
grid on;
xlabel('x Position (m)');
ylabel('y Position (m)');
title(linkTitle);

% Leave room around the link for arrows and labels
axis padded;

limits = axis;
xOrigin = limits(1) + 0.12*(limits(2)-limits(1));
yOrigin = limits(3) + 0.12*(limits(4)-limits(3));
axisLength = 0.15*min(limits(2)-limits(1),limits(4)-limits(3));

quiver(xOrigin,yOrigin,axisLength,0,0,'k','LineWidth',1.5,'MaxHeadSize',0.5);
quiver(xOrigin,yOrigin,0,axisLength,0,'k','LineWidth',1.5,'MaxHeadSize',0.5);

text(xOrigin+axisLength,yOrigin,' +x','FontWeight','bold');
text(xOrigin,yOrigin+axisLength,' +y','FontWeight','bold');

end
