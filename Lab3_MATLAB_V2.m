clc;
clear;

% Lumping units: lengths in mm, masses in g, inertias in g*mm^2.
% Calibrated spring stiffnesses are in N/mm.
% A future dynamic model must convert these quantities to consistent units.
% define variables pressure plate links
L0 = 6; %tenths of an inch this is L1 in the schematic referring to the short distance between lever end and mass
L1 = 40; %tenths of an inch this is L4, 5 in the schematic
L2 = 80; %tenths of an inch this is L2 in the schematic
L3 = 45; %mm
person = 100000;
minput = person/2; %mass assigned to one side; each E/F receives person/4
% scale topout one rev is 300 lbs therefore we get 1.2 degrees per pound.
% Exactly 220 lb corresponds to 264 degrees; 100 kg to about 264.555 degrees.

%pressure plate mechanism variables

L6 = 10; %mm
L7 = 10; %mm
% rack and pinion gearing is 10 tooth pinion 33 tooth rack
% rack has an overall length of 53tenths of an inch with -1tenths of an inch for where the spring hooks
% up. 

% Component masses are omitted provisionally. Zeroing can remove their static
% weight offset, but their inertia still affects a future transient model.

%equation vomit below

%before pressure plate
Jeq1R = (minput/2)*L0^2;
Jeq2R = (minput/2)*L1^2;
Meq3T = Jeq2R/L3^2;
Jeq4R = Meq3T*L1^2;
Meq5T = (Jeq1R + Jeq4R)/L2^2;
Meq6T = 2*Meq5T;

%after pressure plate
M3 = 0; %its zero for now because we do not know the platform mass
Meq7T = Meq6T + M3;
Jeq8R = Meq7T*L6^2;
MR = 0; %it's zero for now, we do not know rack mass
Meq9T = (Jeq8R/(L7^2)) + MR;

%Pinion stuff
Rp = 2; %mm
Jp = 0.5*(1*1^2); %MMOI straightforward, it's a cylinder 
Jeq = Meq9T*Rp^2+Jp; 

%% Spring calibration from expected behavior
% k3: main compression spring at B. k4: weak rack return spring.
% One calibration point identifies kRack = k3*(L6/L7)^2 + k4 ONLY.
% An assumed ratio is needed to split this into two spring constants.
% These are fitted model parameters, not independently measured stiffnesses.
calibrationLoad_lb = 220;
calibrationAngle_deg = 264;
returnSpringRatio = 0.001; %ASSUMPTION: k4/k3 = 0.1%; adjustable, not measured
assert(returnSpringRatio > 0 && returnSpringRatio < 1, ...
    'Return spring stiffness ratio must be between zero and one.');
assert(all([L2,L3,L6,L7,Rp,calibrationLoad_lb,calibrationAngle_deg] > 0), ...
    'Calibration load, angle and transmission lengths must be positive.');
g = 9.80665; %m/s^2
leverRatio = L6/L7; %yB/xRack, under the effective-arm approximation
% Four equal support loads; the two sides give this total force ratio FB/W.
forceRatioB = 0.5*(L0/L2 + (L1/L3)*(L1/L2));
Wcalibration = calibrationLoad_lb*0.45359237*g; %N
FBcalibration = forceRatioB*Wcalibration; %N
xCalibration = Rp*deg2rad(calibrationAngle_deg); %mm
% FB*a = (k3*a^2+k4)*x; substitute k4 = returnSpringRatio*k3.
kRackCalibration = FBcalibration*leverRatio/xCalibration; %N/mm
k3 = kRackCalibration/(leverRatio^2 + returnSpringRatio); %N/mm
k4 = returnSpringRatio*k3; %N/mm
degreesPerPound = calibrationAngle_deg/calibrationLoad_lb;
fprintf('\nSPRING CALIBRATION: %.3f lb -> %.3f deg\n', ...
    calibrationLoad_lb,calibrationAngle_deg);
fprintf('Assumed k4/k3: %.6f (not determined by the calibration point)\n', ...
    returnSpringRatio);
fprintf('k3 = %.6f N/mm; k4 = %.6f N/mm; effective rack stiffness = %.6f N/mm\n', ...
    k3,k4,kRackCalibration);
fprintf('Matching the target is calibration, not independent validation.\n');

%% Static pinion displacement (increment from the unloaded, zeroed position)
% Ideal, frictionless, small-motion model; no damping or transient calculation.
% The listed lengths are treated as effective perpendicular moment arms.
% C is treated as an ideal force-transfer connection between the two levers.
% L1 represents BOTH handwritten L4 (D to F) and L5 (A to C).
% Four equal support loads are assumed. Spring k3 acts at B, k4 at the rack.
% Spring preload and dead weight are absorbed into the unloaded zero setting.
% Gravity forces use ordinary lever ratios, NOT equivalent mass times gravity.
W = (person/1000)*g; %N; convert grams to kilograms
FE = W/4; %N, downward load at E on each side
FF = W/4; %N, downward load at F on each side
FC = FF*L1/L3; %N; moment balance about D
FB_side = (FE*L0 + FC*L1)/L2; %N; moment balance about A
FB = 2*FB_side; %N; both sides load the pressure mechanism

% yB = (L6/L7)*xRack from the intermediate lever's small-motion kinematics.
kRack = k3*leverRatio^2 + k4; %N/mm, stiffness referred to rack
FRackApplied = FB*leverRatio; %N, generalized load referred to rack
xRack = FRackApplied/kRack; %mm
yB = leverRatio*xRack; %mm
thetaPinion_rad = xRack/Rp; %radians; rack travel / pitch radius
thetaPinion_deg = rad2deg(thetaPinion_rad); %total rotation, not modulo 360
input_lb = (person/1000)/0.45359237;
targetAngle_deg = input_lb*degreesPerPound;
indicated_lb = thetaPinion_deg/degreesPerPound;
Fspring3 = k3*yB; %N, incremental spring force
Fspring4 = k4*xRack; %N
FBlever = FB-Fspring3; %N, load carried by intermediate lever

fprintf('\nSTATIC LINEAR MODEL (relative to unloaded zero)\n');
fprintf('Applied mass: %.3f lb; force at B: %.3f N\n', input_lb, FB);
fprintf('Rack travel: %.3f mm; displacement at B: %.3f mm\n', xRack, yB);
fprintf('Pinion rotation: %.3f deg (%.3f rad)\n', thetaPinion_deg, thetaPinion_rad);
fprintf('Calibration target for this mass: %.3f deg\n', targetAngle_deg);
fprintf('Ideal dial reading: %.3f lb\n', indicated_lb);
if max(abs([yB/L6, xRack/L7])) > 0.1
    warning('Scale:SmallMotionExceeded', ...
        ['Predicted travel is large relative to the intermediate lever arms. ', ...
         'Calibration does not validate the small-motion geometry. ', ...
         'Check lever geometry and travel limits before interpreting physical motion.']);
end

%% Free-body diagrams (schematic, not to scale; arrow lengths not force-scaled)
% Individual bodies are isolated to show equal/opposite connection forces.
% Reaction arrows are assumed positive directions; negative results reverse them.
% Spring forces are increments from the zeroed state. Bearing friction and
% dial torque are neglected, so static tangential pinion force Ft is zero.
fbdFigure = figure('Name','Bathroom scale free-body diagrams','Color','w', ...
    'Position',[80 80 1450 850]);
tiledlayout(fbdFigure,2,3,'TileSpacing','loose','Padding','loose');

nexttile; hold on;
plot([0 8],[0 0],'k-','LineWidth',3);
plot(0,0,'kv','MarkerFaceColor','k');
fbdLabel(0,-0.35,'A'); fbdLabel(1.2,-0.35,'E');
fbdLabel(4,-0.35,'C'); fbdLabel(8,-0.35,'B');
fbdForce(1.2,1.8,0,-1.8,sprintf('FE = %.2f N',FE));
fbdForce(4,1.8,0,-1.8,sprintf('FC = %.2f N',FC));
fbdForce(8,-2,0,2,sprintf('FB/2 = %.2f N',FB_side));
fbdForce(0,0,0,2,'Ay'); fbdForce(0,0,-1.4,0,'Ax');
fbdLabel(2,-2.6,'AE = L0; AC = L1; AB = L2');
fbdAxes('1. Upper lever AB (one side)');

nexttile; hold on;
plot([0 7],[0 0],'k-','LineWidth',3); plot(0,0,'kv','MarkerFaceColor','k');
fbdLabel(0,-0.35,'D'); fbdLabel(4,-0.35,'F'); fbdLabel(7,-0.35,'C');
fbdForce(4,2,0,-2,sprintf('FF = %.2f N',FF));
fbdForce(7,-2,0,2,sprintf('FC = %.2f N',FC));
fbdForce(0,0,0,2,'Dy'); fbdForce(0,0,-1.4,0,'Dx');
fbdLabel(2,-2.6,'DF = L1; DC = L3');
fbdAxes('2. Lower lever DC (one side)');

nexttile; hold on;
rectangle('Position',[1,-0.3,6,0.6],'LineWidth',2);
fbdForce(4,2,0,-1.7,sprintf('FB = %.2f N',FB));
fbdForce(2,-2,0,1.7,sprintf('k3 yB = %.2f N',Fspring3));
fbdForce(6,-2,0,1.7,sprintf('Lever reaction = %.2f N',FBlever));
fbdLabel(3,-2.8,'Pressure plate B');
fbdAxes('3. Pressure plate (both sides combined)');

nexttile; hold on;
plot([0 5 5], [1 1 -1.5],'k-','LineWidth',3);
plot(5,1,'kv','MarkerFaceColor','k');
fbdLabel(0,0.6,'B'); fbdLabel(5,1.35,'P'); fbdLabel(5,-1.9,'S');
fbdForce(0,3,0,-2,sprintf('FB - k3 yB = %.2f N',FBlever));
fbdForce(5,-1.5,-2,0,sprintf('Rack reaction = %.3f N',Fspring4));
fbdForce(5,1,1.8,0,'Py'); fbdForce(5,1,0,1.8,'Pz');
fbdLabel(0,-2.8,'PB = L6; PS = L7');
fbdAxes('4. Intermediate lever B-P-S');

nexttile; hold on;
rectangle('Position',[1,-0.25,6,0.5],'LineWidth',2);
fbdForce(-1,0,2,0,sprintf('Lever force = %.3f N',Fspring4));
fbdForce(7,0,-2,0,sprintf('k4 xRack = %.3f N',Fspring4));
fbdLabel(2,1.8,'Ft = 0 at static equilibrium');
fbdLabel(2,-1.7,'Guide reactions omitted (ideal horizontal motion)');
fbdLabel(2,-2.6,'Rack; positive travel to the right');
fbdAxes('5. Rack (incremental horizontal forces)');

nexttile; hold on;
phi = linspace(0,2*pi,200);
plot(4+1.5*cos(phi),0.5+1.5*sin(phi),'k-','LineWidth',2);
plot(4,0.5,'ko','MarkerFaceColor','k');
plot([4 5.5],[0.5 0.5],'k--'); fbdLabel(4.3,0.85,'Rp');
fbdLabel(1,-1.5,'Ft = 0; bearing/dial resisting torque = 0');
fbdLabel(1,-2.2,'theta = xRack / Rp (radians)');
fbdLabel(1,-2.9,sprintf('Linear estimate: %.1f deg',thetaPinion_deg));
fbdAxes('6. Pinion (ideal unloaded dial)');
sgtitle({'Bathroom scale: static free-body diagrams', ...
    'Assumed effective moment arms; incremental forces; schematic geometry'});

fprintf('Ideal dial reading: %.3f lb\n', indicated_lb);
if max(abs([yB/L6, xRack/L7])) > 0.1
    warning('Scale:SmallMotionExceeded', ...
        ['Predicted travel is large relative to the intermediate lever arms. ', ...
        'Calibration does not validate the small-motion geometry. ', ...
        'Check lever geometry and travel limits before interpreting physical motion.']);
end

%% Report plots: static model response over applied weight range
% Uses the same calibrated static model already defined above.
% These plots are for report presentation and do not change the model.

weights_lb = 0:10:300;

theta_deg_plot = zeros(size(weights_lb));
xRack_plot = zeros(size(weights_lb));

FE_plot = zeros(size(weights_lb));
FC_plot = zeros(size(weights_lb));
FB_plot = zeros(size(weights_lb));
FRack_plot = zeros(size(weights_lb));

for i = 1:length(weights_lb)

    % Applied weight
    mass_kg_i = weights_lb(i) * 0.45359237;
    W_i = mass_kg_i * g;

    % Top-view load distribution
    FE_i = W_i/4;
    FF_i = W_i/4;

    FC_i = FF_i * L1/L3;

    FB_side_i = (FE_i*L0 + FC_i*L1)/L2;
    FB_i = 2*FB_side_i;

    % Lever and rack
    FRack_i = FB_i * leverRatio;

    xRack_i = FRack_i / kRack;

    % Pinion / dial rotation
    thetaPinion_rad_i = xRack_i / Rp;
    thetaPinion_deg_i = rad2deg(thetaPinion_rad_i);

    % Store values
    theta_deg_plot(i) = thetaPinion_deg_i;
    xRack_plot(i) = xRack_i;

    FE_plot(i) = FE_i;
    FC_plot(i) = FC_i;
    FB_plot(i) = FB_i;
    FRack_plot(i) = FRack_i;
end


%% Figure: Predicted dial rotation vs applied weight
figure('Color','w');

plot(weights_lb,theta_deg_plot,'LineWidth',1.8);
hold on;

plot(calibrationLoad_lb,calibrationAngle_deg,'o', ...
    'MarkerSize',8,'LineWidth',1.5);

xlabel('Applied Weight (lb)');
ylabel('Pinion / Dial Rotation (deg)');
title('Predicted Dial Rotation vs. Applied Weight');

legend('Static Model', ...
       'Calibration Point', ...
       'Location','northwest');

grid on;


%% Figure: Predicted rack displacement vs applied weight
figure('Color','w');

plot(weights_lb,xRack_plot,'LineWidth',1.8);

xlabel('Applied Weight (lb)');
ylabel('Rack Displacement');
title('Predicted Rack Displacement vs. Applied Weight');

grid on;


%% Figure: Internal force transmission vs applied weight
figure('Color','w');

plot(weights_lb,FE_plot,'LineWidth',1.8);
hold on;

plot(weights_lb,FC_plot,'--','LineWidth',1.8);

plot(weights_lb,FB_plot,'-.','LineWidth',1.8);

xlabel('Applied Weight (lb)');
ylabel('Force (N)');
title('Internal Force Transmission vs. Applied Weight');

legend('Corner Load F_E', ...
       'Transferred Force F_C', ...
       'Platform Force F_B', ...
       'Location','northwest');

grid on;

%% Local plotting helpers
function fbdForce(x,y,dx,dy,label)
quiver(x,y,dx,dy,0,'Color',[0.75 0.12 0.10], ...
    'LineWidth',1.6,'MaxHeadSize',0.35);
label = strrep(label,' = ',sprintf('\n'));
if dy == 0
    text(x+dx/2,y+0.35,label,'FontSize',9,'Interpreter','none', ...
        'HorizontalAlignment','center','VerticalAlignment','bottom');
else
    text(x+0.15,y+dy/2+0.15,label,'FontSize',9,'Interpreter','none', ...
        'VerticalAlignment','bottom');
end
end

function fbdLabel(x,y,label)
text(x,y,label,'FontSize',9,'Interpreter','none');
end

function fbdAxes(titleText)
axis equal; xlim([-2 10]); ylim([-3.5 4]); axis off;
title(titleText,'FontSize',11,'Interpreter','none');
end
