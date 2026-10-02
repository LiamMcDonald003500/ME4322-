clc;
clear;
format short e;
sympref('FloatingPointOutput',true);

%% X-ACTO 1780 mechanical model (requires Symbolic Math Toolbox)

%% Known gearing and observations
teeth = [10 32 12 48 18 56];
N = sym(32)/10*sym(48)/12*sym(56)/18;
a = [sym(1); -sym(10)/32; sym(10)/32*sym(12)/48; -1/N];
n_motor_free = 11000; % rpm
n_cutter_observed = 180; % rpm
omega_motor_free = n_motor_free*2*pi/60;
omega_cutter_observed = n_cutter_observed*2*pi/60;
n_cutter_predicted = n_motor_free/double(N); % rpm
n_input_inferred = n_cutter_observed*double(N);

%% ABS gear mass/inertia estimates from measured geometry
rho_ABS = 1040; % kg/m^3
photo_diameter_cm = [2.35;1.25;0.50;1.80;0.70];
photo_gear_teeth = [56;32;12;48;18];
photo_diameter_m = photo_diameter_cm*0.01;
photo10_measurement_m = 4e-3; % m
syms gear_Do gear_Di gear_thickness positive
gear_volume_annulus = pi/4*(gear_Do^2-gear_Di^2)*gear_thickness;
gear_mass_annulus = rho_ABS*gear_volume_annulus; % kg
gear_weight_annulus = gear_mass_annulus*9.80665; % N
gear_inertia_annulus = gear_mass_annulus*(gear_Do^2+gear_Di^2)/8;
d10 = 5e-3; % m
syms b10 b32 b12 b48 b18 b56 positive
assumed_gear_thickness_m = 3e-3;
gear_OD = [d10 photo_diameter_m(2) photo_diameter_m(3) ...
    photo_diameter_m(4) photo_diameter_m(5) photo_diameter_m(1)];
gear_bore = [b10 b32 b12 b48 b18 b56];
gear_thicknesses = repmat(assumed_gear_thickness_m,1,6);
gear_volumes = pi/4*(gear_OD.^2-gear_bore.^2).*gear_thicknesses;
gear_masses = rho_ABS*gear_volumes;
gear_weights = gear_masses*9.80665;
gear_inertias = gear_masses.*(gear_OD.^2+gear_bore.^2)/8;
J0_ABS_envelope = gear_inertias(1);
J1_ABS_envelope = sum(gear_inertias(2:3));
J2_ABS_envelope = sum(gear_inertias(4:5));
J3_ABS_gear_only = gear_inertias(6);
known_D = double(gear_OD(:));
known_mass_kg = rho_ABS*pi/4*known_D.^2*assumed_gear_thickness_m;
known_J = known_mass_kg.*known_D.^2/8;
motor_proxy.name = "Pololu 1117 generic 130-size brushed DC motor";
motor_proxy.source = "https://www.pololu.com/product/1117/specs";
motor_proxy.voltage_V = 6; % V
motor_proxy.no_load_rpm = 11500;
motor_proxy.shaft_diameter_m = 2e-3;
motor_proxy.case_dimensions_m = [25 15 20]*1e-3;
motor_proxy.total_mass_kg = 18e-3; % kg
motor_proxy.no_load_current_A = 0.070;
motor_proxy.stall_current_A = 0.800;
motor_proxy.speed_difference_percent = ...
    100*(motor_proxy.no_load_rpm-n_motor_free)/n_motor_free;

%% Unknown dynamic parameters: retain symbols until geometry is confirmed
syms Jr J0 J1 J2 J3 Ks positive
syms Ds Dr D0 D1 D2 D3 nonnegative
syms Tmotor Tblade real
Jg = simplify(J0+a(2)^2*J1+a(3)^2*J2+a(4)^2*J3);
syms J0_extra J1_extra J2_extra J3_extra nonnegative
assembly_J_geometry = [J0_ABS_envelope+J0_extra; ...
    J1_ABS_envelope+J1_extra;J2_ABS_envelope+J2_extra; ...
    J3_ABS_gear_only+J3_extra];
Jg_geometry = simplify(subs(Jg,[J0 J1 J2 J3],assembly_J_geometry.'));
Dg = simplify(D0+a(2)^2*D1+a(3)^2*D2+a(4)^2*D3);
syms J0_remaining J1_remaining J2_remaining J3_remaining real
disk_lump_inertias = [known_J(1);sum(known_J(2:3)); ...
    sum(known_J(4:5));known_J(6)];
disk_reflected_inertia = sum(double(a).^2.*disk_lump_inertias);
remaining_inertias = [J0_remaining;J1_remaining;J2_remaining;J3_remaining];
reflected_remaining = simplify(sum(a.^2.*remaining_inertias));
assembly_J_disk_with_remaining = disk_lump_inertias+remaining_inertias;

%% Geometry formulas for later inertia and stiffness estimates
syms mass radius radius_inner radius_outer offset Jcm real
syms Gshaft shaft_diameter shaft_length positive
J_solid = mass*radius^2/2;
J_annular = mass*(radius_outer^2+radius_inner^2)/2;
J_offset = Jcm+mass*offset^2;
Jpolar = pi*shaft_diameter^4/32; % m^4
Ks_geometry = Gshaft*Jpolar/shaft_length;

syms module1 module2 module3 positive
gear_modules = [module1 module1 module2 module2 module3 module3];
pitch_diameters = gear_modules.*teeth;
outside_diameters = gear_modules.*(teeth+2);
pitch_radii = pitch_diameters/2;
circular_pitches = pi*gear_modules;
mesh_center_distances = [module1*(10+32)/2; ...
    module2*(12+48)/2;module3*(18+56)/2];

%% Equations of motion
syms theta_r q omega_r omega_g alpha_r alpha_g real
angles = [theta_r;q];
velocities = [omega_r;omega_g];
accelerations = [alpha_r;alpha_g];
shaft_twist = theta_r-q;
shaft_twist_rate = omega_r-omega_g;
shaft_torque = Ks*shaft_twist+Ds*shaft_twist_rate;
cutting_torque_input = Tblade/N;
M = diag([Jr Jg]);
D = [Dr+Ds -Ds;-Ds Dg+Ds];
K = [Ks -Ks;-Ks Ks];
input_map = [sym(1) sym(0);sym(0) -1/N];
u = [Tmotor;Tblade];
force = input_map*u;
eom = M*accelerations+D*velocities+K*angles == force;

%% First-order state equations: x=[theta_r;q;omega_r;omega_g]
x = [angles;velocities];
A = [sym(zeros(2)) sym(eye(2));-(M\K) -(M\D)];
Bu = [sym(zeros(2));M\input_map];
xdot = simplify(A*x+Bu*u);

%% Outputs and gear kinematics
Cout = [0 0 1 0;0 0 0 -1/N;1 -1 0 0;Ks -Ks Ds -Ds];
Dout = sym(zeros(4,2));
y = simplify(Cout*x+Dout*u);
output_names = ["Motor speed (rad/s)";"Cutter speed (rad/s)";"Shaft twist (rad)";"Shaft torque (N*m)"];
gear_angles = a*q;
gear_speeds = a*omega_g;
gear_accelerations = a*alpha_g;
theta_out = -q/N;
omega_out = -omega_g/N;
motor_rpm = omega_r*60/(2*pi);
cutter_rpm = omega_out*60/(2*pi); % rpm

%% Initial conditions: assumed startup from rest, zero initial twist
x0 = zeros(4,1);

%% Symbolic Laplace solutions and transfer matrix
syms s Tmotor_s Tblade_s real
syms theta_r0 q0 omega_r0 omega_g0 real
dynamic_matrix = M*s^2+D*s+K;
theta_s_zeroIC = simplify(dynamic_matrix\(input_map*[Tmotor_s;Tblade_s]));
theta0 = [theta_r0;q0];
omega0 = [omega_r0;omega_g0];
theta_s_general = simplify(dynamic_matrix\(input_map*[Tmotor_s;Tblade_s]+M*(s*theta0+omega0)+D*theta0));
angle_transfer = simplify(dynamic_matrix\input_map);
H = simplify([s 0;0 -s/N;1 -1;Ks+Ds*s -Ks-Ds*s]*angle_transfer);

%% Symbolic constant-load steady rotation and torsional mode
steady_speed_input = simplify((Tmotor-Tblade/N)/(Dr+Dg));
steady_speed_output = -steady_speed_input/N;
steady_shaft_torque = simplify(Dg*steady_speed_input+Tblade/N);
steady_shaft_twist = simplify(steady_shaft_torque/Ks);
steady_no_pencil_speed = simplify(subs(steady_speed_input,Tblade,0));
torsional_natural_frequency = sqrt(Ks*(1/Jr+1/Jg)); % rad/s
torsional_natural_frequency_Hz = torsional_natural_frequency/(2*pi);
characteristic_polynomial = collect(expand(det(dynamic_matrix)),s);


%% Energy balance verification
energy = Jr*omega_r^2/2+Jg*omega_g^2/2+Ks*shaft_twist^2/2;
energy_rate = Tmotor*omega_r-Dr*omega_r^2-Dg*omega_g^2 ...
    -Ds*shaft_twist_rate^2-Tblade*omega_g/N;
energy_balance_residual = simplify(jacobian(energy,x)*xdot-energy_rate);
assert(isAlways(energy_balance_residual==0),'Energy balance failed.');

%% Propagate confirmed geometry into the dynamic model
inertia_symbols = [J0 J1 J2 J3];
inertia_geometry = assembly_J_geometry.';
M_geometry = subs(M,inertia_symbols,inertia_geometry);
A_geometry = subs(A,inertia_symbols,inertia_geometry);
A_disk_with_remaining = subs(A,inertia_symbols,assembly_J_disk_with_remaining.');
Bu_geometry = subs(Bu,inertia_symbols,inertia_geometry);
eom_geometry = subs(eom,inertia_symbols,inertia_geometry);
xdot_geometry = subs(xdot,inertia_symbols,inertia_geometry);
H_geometry = subs(H,inertia_symbols,inertia_geometry);
theta_s_geometry = subs(theta_s_general,inertia_symbols,inertia_geometry);
energy_geometry = subs(energy,inertia_symbols,inertia_geometry);
natural_frequency_geometry = subs(torsional_natural_frequency, ...
    inertia_symbols,inertia_geometry);
characteristic_polynomial_geometry = subs(characteristic_polynomial, ...
    inertia_symbols,inertia_geometry);
syms J0_correction J1_correction J2_correction J3_correction real
assembly_J_corrected = assembly_J_geometry + ...
    [J0_correction;J1_correction;J2_correction;J3_correction];
Jg_corrected = subs(Jg,inertia_symbols,assembly_J_corrected.');


%% Required analysis results only
%% Pencil inserted: single-blade cutting load and blade inertia
blade_length_m = 18.5e-3;
blade_width_m = 8e-3;
blade_thickness_m = 0.35e-3;
syms rho_blade positive
syms blade_offset nonnegative
syms axis_L axis_W axis_H real
blade_volume_m3 = blade_length_m*blade_width_m*blade_thickness_m;
blade_mass = rho_blade*blade_volume_m3;
blade_J_centroid = blade_mass/12*( ...
    axis_L^2*(blade_width_m^2+blade_thickness_m^2) + ...
    axis_W^2*(blade_length_m^2+blade_thickness_m^2) + ...
    axis_H^2*(blade_length_m^2+blade_width_m^2));
blade_J = blade_J_centroid+blade_mass*blade_offset^2;
cutter_length_m = 2.5e-2;
cutter_width_m = 1e-2;
cutter_height_m = 1.5e-2;
rho_PS = 1040; % kg/m^3
cutter_body_mass = rho_PS*cutter_length_m*cutter_width_m*cutter_height_m;
cutter_body_weight = cutter_body_mass*9.80665;
cutter_body_J = cutter_body_mass*(cutter_width_m^2+cutter_height_m^2)/12;
syms Joutput_other nonnegative
syms Jbody_correction real
cut_inertias = subs(assembly_J_corrected.',J3_extra, ...
    cutter_body_J+Jbody_correction+Joutput_other+blade_J);
Jg_cut = subs(Jg,inertia_symbols,cut_inertias);
syms Tcut Dc nonnegative
cutting_load = Tcut+Dc*omega_g/N; % N*m
cut_eom = subs(eom,inertia_symbols,cut_inertias);
cut_eom = subs(cut_eom,Tblade,cutting_load);
A_cut = subs(A,inertia_symbols,cut_inertias);
A_cut(4,4) = A_cut(4,4)-Dc/(N^2*Jg_cut);
Bu_cut = subs(Bu,inertia_symbols,cut_inertias);
xdot_cut = A_cut*x+Bu_cut*[Tmotor;Tcut];
H_cut = subs(subs(H,inertia_symbols,cut_inertias),D0,D0+Dc/N^2);
steady_speed_cut = (Tmotor-Tcut/N)/(Dr+Dg+Dc/N^2);
steady_cutting_load = Tcut+Dc*steady_speed_cut/N;
steady_torque_cut = Dg*steady_speed_cut+steady_cutting_load/N;
steady_twist_cut = steady_torque_cut/Ks;
steady_outputs_cut = [steady_speed_cut;-steady_speed_cut/N; ...
    steady_torque_cut;steady_twist_cut];
syms Tcut_s real
theta_s_cut = subs(subs(subs(theta_s_general,inertia_symbols,cut_inertias), ...
    D0,D0+Dc/N^2),Tblade_s,Tcut_s);
natural_frequency_cut = sqrt(Ks*(1/Jr+1/Jg_cut));
energy_cut = subs(energy,inertia_symbols,cut_inertias);
energy_rate_cut = Tmotor*omega_r-Dr*omega_r^2-Dg*omega_g^2 ...
    -Ds*shaft_twist_rate^2-Tcut*omega_g/N-Dc*omega_g^2/N^2;
cut_energy_residual = simplify(jacobian(energy_cut,x)*xdot_cut-energy_rate_cut);
assert(isAlways(cut_energy_residual==0),'Pencil-inserted energy balance failed.');

%% Report-ready lumping derivation and final differential equations
syms J10 J32 J12 J48 J18 J56 Jshaft0 Jshaft1 Jshaft2 Jshaft3 nonnegative
component_lumps = [J10+Jshaft0; J32+J12+Jshaft1; ...
    J48+J18+Jshaft2; J56+Jshaft3+cutter_body_J+Jbody_correction+blade_J];
kinetic_energy_gears = sum([J0;J1;J2;J3].*(a*omega_g).^2)/2;
reflected_damping_terms = a.^2.*[D0;D1;D2;D3];
reflected_inertia_terms = a.^2.*cut_inertias.';
input_cutting_damping = Dc/N^2;
M_cut = subs(M,inertia_symbols,cut_inertias);
D_cut = D; D_cut(2,2) = D_cut(2,2)+input_cutting_damping;
K_cut = K;
force_cut = input_map*[Tmotor;Tcut];
cut_dynamic_residual = simplify(M_cut*xdot_cut(3:4)+ ...
    D_cut*velocities+K_cut*angles-force_cut);
assert(all(isAlways(cut_dynamic_residual==0)), ...
    'Final differential equation disagrees with the state model.');
assert(isAlways(simplify(kinetic_energy_gears-Jg*omega_g^2/2)==0), ...
    'Squared-speed inertia reflection failed.');

syms Jg_eff positive
syms Dg_eff nonnegative
syms t real
syms theta_r_t(t) q_t(t) Tmotor_t(t) Tcut_t(t)
rotor_ode = Jr*diff(theta_r_t(t),t,2)+(Dr+Ds)*diff(theta_r_t(t),t) ...
    -Ds*diff(q_t(t),t)+Ks*(theta_r_t(t)-q_t(t)) == Tmotor_t(t);
gearbox_ode = Jg_eff*diff(q_t(t),t,2)-Ds*diff(theta_r_t(t),t) ...
    +(Dg_eff+Ds+Dc/N^2)*diff(q_t(t),t) ...
    +Ks*(q_t(t)-theta_r_t(t)) == -Tcut_t(t)/N;
final_differential_equations = [rotor_ode;gearbox_ode];
M_report = diag([Jr Jg_eff]);
D_report = [Dr+Ds -Ds;-Ds Dg_eff+Ds+Dc/N^2];
A_report = [sym(zeros(2)) sym(eye(2)); ...
    -(M_report\K) -(M_report\D_report)];
Bu_report = [sym(zeros(2));M_report\input_map];
dynamic_report = M_report*s^2+D_report*s+K;
H_report = simplify([s 0;0 -s/N;1 -1;Ks+Ds*s -Ks-Ds*s] ...
    *(dynamic_report\input_map));
theta_s_report = simplify(dynamic_report\(input_map*[Tmotor_s;Tcut_s] ...
    +M_report*(s*theta0+omega0)+D_report*theta0));
omega_steady_report = (Tmotor-Tcut/N)/(Dr+Dg_eff+Dc/N^2);
torque_steady_report = (Dg_eff+Dc/N^2)*omega_steady_report+Tcut/N;
steady_report = [omega_steady_report;-omega_steady_report/N; ...
    torque_steady_report;torque_steady_report/Ks];
frequency_report = sqrt(Ks*(1/Jr+1/Jg_eff));

%% Required parameter table: values, units, status, measurement method
parameter_names = ["Jr";"J0";"J1";"J2";"J3";"Jg_eff"; ...
    "Ks";"Ds";"Dr";"D0";"D1";"D2";"D3";"Dg_eff";"Dc"; ...
    "Tmotor";"Tcut";"N";"cutter_body_mass";"cutter_body_J";"blade_J"];
parameter_expressions = [{Jr};num2cell(cut_inertias.'); ...
    {Jg_cut;Ks;Ds;Dr;D0;D1;D2;D3;Dg;Dc;Tmotor;Tcut;N; ...
     cutter_body_mass;cutter_body_J;blade_J}];
parameter_units = [repmat("kg*m^2",6,1);"N*m/rad"; ...
    repmat("N*m*s/rad",8,1);"N*m";"N*m";"1";"kg";"kg*m^2";"kg*m^2"];
parameter_status = ["Unknown rotor inertia"; ...
    repmat("Geometry estimate with symbolic remaining parts/corrections",4,1); ...
    "Derived from shaft inertias and squared speed ratios"; ...
    repmat("Unknown: mechanical measurement required",9,1); ...
    "Prescribed motor-torque input";"Unknown pencil-cutting load"; ...
    "Inspected tooth counts";"Assumed solid polystyrene body"; ...
    "Assumed centered longitudinal axis";"Blade material and mounting unresolved"];
parameter_method = ["Torsional pendulum or rotor geometry/material model"; ...
    repmat("Measure bores/hubs/shafts; subtract tooth gaps and sum inertias",4,1); ...
    "Jg_eff = sum(a_i^2*J_i), with cutter included in J3"; ...
    "Ks = G*pi*d^4/(32*L), or measured torque/twist"; ...
    "Fit shaft relative-motion decay or torque versus relative speed"; ...
    repmat("Isolate rotor/shaft; fit unloaded coast-down with known inertia",5,1); ...
    "Dg_eff = sum(a_i^2*D_i)"; ...
    "Fit added cutting torque versus cutter speed at fixed pencil feed"; ...
    "Mechanical torque fixture; speed alone does not establish torque"; ...
    "Measure cutting reaction torque at controlled pencil feed"; ...
    "(32/10)*(48/12)*(56/18)"; ...
    "rho_PS*length*width*height; weigh or correct for holes when available"; ...
    "mass*(width^2+height^2)/12; correct for removed material"; ...
    "Measure density, axis orientation and offset; J = J_axis,centroid + m*offset^2"];
parameter_values = strings(numel(parameter_names),1);
for row = 1:numel(parameter_names)
    expression = sym(parameter_expressions{row});
    if isempty(symvar(expression))
        parameter_values(row) = string(sprintf('%.6e',double(expression)));
    else
        parameter_values(row) = string(char(vpa(expression,5)));
    end
end
parameter_results = table(parameter_names,parameter_values,parameter_units, ...
    parameter_status,parameter_method);
analysis_output_dir = fullfile(fileparts(mfilename('fullpath')),'analysis_output');
if ~exist(analysis_output_dir,'dir'), mkdir(analysis_output_dir); end
writetable(parameter_results,fullfile(analysis_output_dir,'model_parameters.csv'));

%% Required analysis results only
fprintf('GEAR TRANSMISSION AND SPEED CHECKS\n');
fprintf('Reduction: %.4e:1\n',double(N));
fprintf('No-pencil conditional output at 11000 rpm input: %.4e rpm\n',n_cutter_predicted);
fprintf('Input inferred from no-pencil 180 rpm output: %.4e rpm\n\n',n_input_inferred);
fprintf('GEAR-ONLY MMOI ESTIMATES (kg*m^2; ABS disks, 3 mm thick)\n');
for k = 1:4
    fprintf('J%d gear contribution: %.6e\n',k-1,disk_lump_inertias(k));
end
fprintf('Reflected gear-only Jg: %.6e\n\n',disk_reflected_inertia);
disp('MODEL ASSUMPTIONS AND JUSTIFICATION');
disp('Fixed housing and held pencil: model rotational drivetrain motion only.');
disp('Rigid compound gears and negligible backlash: one speed per shaft.');
disp('Flexible motor shaft: retains relative twist and transmitted torque.');
disp('Viscous bearing/lubrication loss: first approximation, identified mechanically.');
disp('Cutting resistance = Tcut + Dc*abs(cutter speed), for forward cutting.');
disp('ABS tooth-tip envelopes and solid PS body estimate geometry; corrections remain symbolic.');
disp('No independent translating mass: mass contributes through rotational MMOI.');

disp('LUMPING STEP 1: INDIVIDUAL PARTS -> RIGID SHAFT ASSEMBLIES [J0; J1; J2; J3]');
disp(vpa(component_lumps,5));
disp('LUMPING STEP 2: SIGNED SHAFT SPEEDS [omega0; omega1; omega2; omega3] = a*omega_g');
for k = 1:4
    fprintf('a%d = %s\n',k-1,char(a(k)));
end
disp('LUMPING STEP 3: kinetic energy = 1/2*sum(J_i*(a_i*omega_g)^2) = 1/2*Jg*omega_g^2');
disp('Each input-reflected inertia contribution, with supplied geometry:');
disp(vpa(reflected_inertia_terms,5));
disp('Jg_eff = sum of these contributions:');
disp(vpa(Jg_cut,5));
disp('LUMPING STEP 4: bearing dissipation uses squared speed ratios');
disp('Input-reflected [D0; a1^2 D1; a2^2 D2; a3^2 D3]:');
disp(reflected_damping_terms);
disp('Dg_eff = sum of these terms:');
disp(Dg);
disp('LUMPING STEP 5: opposing cutter torque reflects as Tcut/N; cutting damping as Dc/N^2');
disp([Tcut/N;input_cutting_damping]);
disp('LUMPED MODEL BLOCK DIAGRAM (rotational mechanical connections)');
disp('Tmotor -> [rotor Jr, Dr; theta_r] -- [shaft Ks || Ds] -- [gearbox Jg_eff; q]');
disp('                                                     | bearing Dg_eff to ground');
disp('                                                     | cutting Dc/N^2 to ground');
disp('                                                     <- opposing load Tcut/N');
disp('Output relation: theta_out = -q/N; omega_out = -qdot/N.');
disp('FINAL PENCIL-INSERTED DIFFERENTIAL EQUATIONS (d/dt notation)');
disp(final_differential_equations);
disp('Jg_eff is the geometry-informed inertia above; Dg_eff is bearing damping above.');
disp('FINAL MATRIX EOM: M_cut*thetaddot + D_cut*thetadot + K_cut*theta = force_cut');
disp('theta = [theta_r; q]. Matrices M, D, K and force follow:');
disp([Jr 0;0 Jg_eff]);
disp([Dr+Ds -Ds;-Ds Dg_eff+Ds+Dc/N^2]);
disp(K_cut);
disp(force_cut);
disp('SHAFT STIFFNESS FROM GEOMETRY (shaft dimensions/material still required)');
disp(Ks_geometry);
disp('REQUIRED PARAMETER TABLE (CSV saved in analysis_output/model_parameters.csv)');
disp(parameter_results(:,{'parameter_names','parameter_units','parameter_status'}));
disp('Full numerical/symbolic values and measurement methods are in the CSV table.');
disp('STATE-SPACE MATRICES: xdot = A*x + Bu*u');
disp(vpa(A_report,5));
disp(vpa(Bu_report,5));
disp('OUTPUT EXPRESSIONS: motor speed, cutter speed, shaft twist, shaft torque');
disp(vpa(y,5));
disp('TRANSFER MATRIX: same output order; inputs are motor torque and Tcut');
disp(vpa(H_report,5));
disp('LAPLACE RESPONSE INCLUDING INITIAL CONDITIONS');
disp(vpa(theta_s_report,5));
disp('STEADY ROTATION: motor/input speed, cutter speed, shaft torque, shaft twist');
disp(vpa(steady_report,5));
disp('UNDAMPED TORSIONAL NATURAL FREQUENCY (rad/s)');
disp(vpa(frequency_report,5));
disp('UNDAMPED TORSIONAL NATURAL FREQUENCY (Hz)');
disp(vpa(frequency_report/(2*pi),5));
disp('ENERGY BALANCE RESIDUAL (must be zero)');
disp(cut_energy_residual);
disp('PENCIL INSERTED: cutter resistance = Tcut + Dc*abs(output speed)');
disp('BLADE MMOI: density, axis direction, offset, and holder inertia remain symbolic');
disp(vpa(blade_J,5));
fprintf('ASSUMED SOLID POLYSTYRENE CUTTER BODY\n');
fprintf('Mass: %.4e kg; weight: %.4e N; longitudinal MMOI: %.4e kg*m^2\n', ...
    cutter_body_mass,cutter_body_weight,cutter_body_J);
fprintf('Startup initial state: [0; 0; 0; 0].\n');
fprintf('theta_r(0)=q(0)=omega_r(0)=omega_g(0)=0; initial shaft twist=0.\n');
fprintf('Time-response plots await remaining parameters and torque histories.\n');

%% Generic motor efficiency comparison (separate from mechanical EOM)
pencil_description = "Ticonderoga pencil with cedar casing (user identified)";
V_proxy = motor_proxy.voltage_V;
I0_proxy = motor_proxy.no_load_current_A;
Ist_proxy = motor_proxy.stall_current_A;
omega0_proxy = motor_proxy.no_load_rpm*2*pi/60;
R_proxy = V_proxy/Ist_proxy;
Kt_proxy = (V_proxy-R_proxy*I0_proxy)/omega0_proxy;
Tstall_net_proxy = Kt_proxy*(Ist_proxy-I0_proxy);
Pshaft_max_proxy = Tstall_net_proxy*omega0_proxy/4;
eta_motor_peak_proxy = (1-sqrt(I0_proxy/Ist_proxy))^2;
n_motor_peak_eta_proxy = motor_proxy.no_load_rpm/(1+sqrt(I0_proxy/Ist_proxy));

syms n_cutter_loaded positive
omega_loaded_proxy = double(N)*n_cutter_loaded*2*pi/60;
I_loaded_proxy = (V_proxy-Kt_proxy*omega_loaded_proxy)/R_proxy;
Tshaft_loaded_proxy = Kt_proxy*(I_loaded_proxy-I0_proxy);
eta_motor_loaded_proxy = Tshaft_loaded_proxy*omega_loaded_proxy/(V_proxy*I_loaded_proxy);
syms eta_gearbox positive
eta_overall_loaded_proxy = eta_motor_loaded_proxy*eta_gearbox;
eta_mechanical_steady = (Tcut/N+Dc*omega_steady_report/N^2)/Tmotor;

eta_mesh_cases = [0.90;0.95;0.98];
eta_gear_cases = eta_mesh_cases.^3;
eta_overall_peak_cases = eta_motor_peak_proxy*eta_gear_cases;
efficiency_scenarios = table(100*eta_mesh_cases,100*eta_gear_cases, ...
    100*eta_overall_peak_cases,'VariableNames', ...
    {'assumed_mesh_efficiency_percent','gearbox_efficiency_percent', ...
    'peak_electrical_to_cutter_efficiency_percent'});
writetable(efficiency_scenarios,fullfile(analysis_output_dir,'efficiency_scenarios.csv'));

gearmotor_reference_stall_Nm = 75*0.00706155183333;
gearmotor_reference_ratio = (36*32*37*24)/(8*14*9*8);
assert(eta_motor_peak_proxy>0 && eta_motor_peak_proxy<1);
assert(abs(double(subs(eta_motor_loaded_proxy,n_cutter_loaded, ...
    n_motor_peak_eta_proxy/double(N)))-eta_motor_peak_proxy)<1e-10);
fprintf('\nEFFICIENCY: GENERIC MOTOR ESTIMATE, NOT MEASURED SHARPENER PERFORMANCE\n');
fprintf('Pencil: %s\n',pencil_description);
fprintf('Estimated R: %.4e ohm; Kt: %.4e N*m/A\n',R_proxy,Kt_proxy);
fprintf('Estimated net stall torque: %.4e N*m; max shaft power: %.4e W\n', ...
    Tstall_net_proxy,Pshaft_max_proxy);
fprintf('Estimated peak motor efficiency: %.2f percent at %.4e motor rpm\n', ...
    100*eta_motor_peak_proxy,n_motor_peak_eta_proxy);
fprintf('Assumed mesh %% | gearbox %% | peak electrical-to-cutter %%\n');
for k=1:numel(eta_mesh_cases)
    fprintf('%.2f | %.2f | %.2f\n',100*eta_mesh_cases(k), ...
        100*eta_gear_cases(k),100*eta_overall_peak_cases(k));
end
disp('Peak cases exclude additional bearing/support losses and battery losses.');
disp('Loaded operating efficiency remains a function of pencil-inserted cutter rpm:');
disp(vpa(eta_overall_loaded_proxy,5));
disp('Mechanical EOM efficiency at positive steady rotation:');
disp(eta_mechanical_steady);
fprintf('75 oz-in = %.4e N*m is a separate %.4e:1 gearmotor OUTPUT stall rating.\n', ...
    gearmotor_reference_stall_Nm,gearmotor_reference_ratio);
