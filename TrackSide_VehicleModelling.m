%% ________________________________________________________________________
%% TRACK SIDE VEHICLE MODELLING SCRIPT
%
% Vehicle Model generation for use in lap-time simulation
%
% Liscensing: GPL V3 Open Source License.
% Information:
%   Nirat Pai
%   BSc Mechanical Engineering, Purdue University
%   Minor in Electrical & Computer Engineering, and Economics
%   LinkedIn: linkedin.com/in/nirat-pai
%   email: niratpai24@gmail.com
%   GitHub: https://github.com/waffless24

%% ________________________________________________________________________
%% CLEARING WORKSPACE
clear
clc

%% ________________________________________________________________________
%% READING FILES

% Set the desired vehicle
vehicle_name = 'Mercedes AMG One';

% Obtaining files from folder
directory_name = append('Vehicles/', vehicle_name, '/');
directory = dir(append('Vehicles/', vehicle_name, '/*.xlsx'));

% Sorting into vehicle components
engine = directory(find(contains({directory.name},'Engine'))).name;
chassis = directory(find(contains({directory.name},'Chassis'))).name;
brakes = directory(find(contains({directory.name},'Brakes'))).name;
tyres = directory(find(contains({directory.name},'Tyres'))).name;

engine_spec = read_specs(directory_name, engine, 'Specs');
torque_curve = read_torque_curve(directory_name, engine, 'Torque Curve');
clear("engine");

brakes = read_specs(directory_name, brakes);
tyres = read_specs(directory_name, tyres);
chassis = read_specs(directory_name, chassis);

% Displaying Vehicle Information
disp('________________________________________________________________________');
fprintf('<strong> %s </strong>\n', vehicle_name);
disp('________________________________________________________________________');
fprintf('<strong> Engine:</strong> %s \n', engine_spec{1,"Value"});
fprintf('<strong> Engine Type:</strong> %s \n', engine_spec{2,"Value"});
fprintf('<strong> Chassis:</strong> %s, %s\n', chassis{1,"Value"}, chassis{2,"Value"});
fprintf('<strong> Brakes:</strong> %s \n', brakes{1,"Value"});
fprintf('<strong> Tyres:</strong> %s \n', tyres{1,"Value"});
disp('________________________________________________________________________');

%% ________________________________________________________________________
%% IMPORTING VARIABLES

% CHASSIS:
M = double(chassis{7, "Value"});
M_dist = double(chassis{8, "Value"});
Cl = double(chassis{9, "Value"});
Cl_DRS = double(chassis{14, "Value"});
Cd = double(chassis{10, "Value"});
Cd_DRS = double(chassis{15, "Value"});
aero_dist = double(chassis{11, "Value"});
A = double(chassis{12, "Value"});

% ENGINE:
power_factor = double(engine_spec{4, "Value"});
thermal_e = double(engine_spec{5, "Value"});
fuel_tank_cap = double(engine_spec{6, "Value"});
fuel_LHV = double(engine_spec{7, "Value"});

if engine_spec{2,"Value"} == "Hybrid"
    battery_cap = engine_spec{8,"Value"};
end

drive_type = engine_spec{9, "Value"};

% TRANSMISSION:
e_driveline = double(engine_spec{10, "Value"});
r_primary = double(engine_spec{11, "Value"});
r_drive = double(engine_spec{12, "Value"});
r_gearbox = double(engine_spec{13:end, "Value"});
nog = length(r_gearbox);
r_final = r_primary * r_drive * r_gearbox;

% BRAKES:
d_disc_front = double(brakes{2, "Value"})/1000;
d_disc_rear = double(brakes{3, "Value"})/1000;
hpad_front = double(brakes{4, "Value"})/1000;
hpad_rear =double(brakes{5, "Value"})/1000;

mu_pad = double(brakes{6, "Value"});

nop_front = double(brakes{7, "Value"});
nop_rear = double(brakes{8, "Value"});

d_cal = double(brakes{9, "Value"})/1000;
d_mast = double(brakes{10, "Value"})/1000;

r_pedal = double(brakes{11, "Value"});

% TYRES:
tyre_rad = (double(tyres{3, "Value"}) * 0.0254) / 2;

% Coefficient of Rolling
mu_r = double(tyres{4, "Value"}); 

% Longitudnal Coefficients
mu_y = double(tyres{5, "Value"}); 
loadR_y = double(tyres{6, "Value"});
mu_y_sens = double(tyres{7, "Value"});

% Lateral Coefficients
mu_x = double(tyres{8, "Value"});
loadR_x = double(tyres{9, "Value"});
mu_x_sens = double(tyres{10, "Value"});

% Tyre stiffness (Torque Resistance)
stiff_F = double(tyres{11, "Value"});
stiff_R = double(tyres{12, "Value"});

%% ________________________________________________________________________
%% POWERTRAIN MODEL

% Engine Torque & Power Curve
engine_RPM = double(torque_curve{:,"RPM"});
engine_torque = double(torque_curve{:,"Torque"});
engine_power = engine_torque .* engine_RPM * (pi/30);

% Finding available torque and vehicle speed for each gear set
gear_speeds = zeros(length(engine_RPM), nog);
gear_torque = zeros(length(engine_RPM), nog);

for i = 1:nog
    gear_speeds(:,i) = (engine_RPM * (pi/30) / r_final(i)) * tyre_rad;
    gear_torque(:,i) = engine_torque * r_final(i) * e_driveline;
end

% Generating Vehicle speed array
v_min =min(gear_speeds, [], "all");
v_max = max(gear_speeds, [], "all");
vehicle_speed = linspace(v_min, v_max , 499).';
gear_tf = zeros(length(vehicle_speed), nog);

for i = 1:nog
    gear_tf(:,i) = interp1(gear_speeds(:,i), gear_torque(:,i)/tyre_rad, vehicle_speed, 'spline', 0);
end

% Selecting Gear based on max available tractive force
[tf_max, gears] = max(gear_tf,[], 2);

% Adding in zeros for better modelling
vehicle_speed  = [0;vehicle_speed];
gears = [gears(1); gears];
gear_tf = [zeros(1,nog);gear_tf];
tf_max = [0; tf_max];

disp('Powertrain Model generated....');

% Gear shift points & Engine RPM drops
engine_RPM_mod = vehicle_speed/tyre_rad .* r_final(gears) * (30/pi);
g_shifts = find(diff(gears));
g_shifts = vehicle_speed(g_shifts);

disp('Gear Shift points genearted....');

%% ________________________________________________________________________
%% BRAKE MODEL

% Dimensions
a_mast = (pi * d_mast^2) / 4;
a_cal_front = (nop_front * pi * d_cal^2) / 4;
a_cal_rear = (nop_rear * pi * d_cal^2) / 4;

% Brake coefficients required for modelling
brc1_front = 0.25 * (tyre_rad / ((d_disc_front - hpad_front)/ 2)) * (1/(mu_pad * a_cal_front));
brc1_rear = 0.25 * (tyre_rad / ((d_disc_rear - hpad_rear)/ 2)) * (1/(mu_pad * a_cal_rear));

brc2 = a_mast / r_pedal;

disp('Brake model genearted....');
%% ________________________________________________________________________
%% STEERING MODEL

%% ________________________________________________________________________
%% FORCE MODEL
% Constants
g = 9.806;
rho = 1.226;

bank = 0;
incl = 0;

% Natural forces
fz_weight = - M * g * cosd(bank) * cosd(incl);
fz_aero = 0.5 * rho * A * Cl * (vehicle_speed.^2);
fz_total = fz_weight + fz_aero;

fy_weight = M * g * sind(incl);
fy_aero = 0.5 * rho * A * Cd * (vehicle_speed.^2);
fy_roll = mu_r * fz_total;

fx_weight = - M * g * sind(bank);

% Tyre Forces
Ny = loadR_y * g;
Nx = loadR_x * g;
fn_front = ((M_dist) * fz_weight + (aero_dist) * fz_aero) / 2;
fn_rear = ((1 - M_dist) * fz_weight + (1 - aero_dist) * fz_aero) / 2;

switch drive_type
    case 'RWD'
        fy_tyres = (mu_y + mu_y_sens * (Ny - abs(fn_rear))) .* abs(fn_rear) * 2;
    case 'FWD'
        fy_tyres = (mu_y + mu_y_sens * (Ny * g - abs(fn_front))) .* abs(fn_front) * 2;
    otherwise %AWD
        fy_rear = (mu_y + mu_y_sens * (Ny * g - abs(fn_rear))) .* abs(fn_rear) * 2;
        fy_front = (mu_y + mu_y_sens * (Ny * g - abs(fn_front))) .* abs(fn_front) * 2;
        fy_tyres = fy_rear + fy_front;
end

fx_rear = (mu_x + mu_x_sens * (Nx * g - abs(fn_rear))) .* abs(fn_rear) * 2;
fx_front = (mu_x + mu_x_sens * (Nx * g - abs(fn_front))) .* abs(fn_front) * 2;
fx_tyres = fx_rear + fx_front;

disp('Force Model generated....');

%% ________________________________________________________________________
%% GGV MAP

% Number of points on the ellipse
ep = 90;

% Initializing GGV matrix
GGV = zeros(length(vehicle_speed), 2 * ep-1, 3);

for i = 1:length(vehicle_speed)
    % Longitudnal resistive acceleration
    ay_resistive = (fy_weight + fy_aero(i) + fy_roll(i)) / M;
    
    % Power Limited acceleration
    ay_engine = tf_max(i) / M;
    ay_engine = ay_engine * ones(ep,1) ;

    % Tyre acceleration
    ay_tyres_acc_max = fy_tyres(i) / M;
    ay_tyres_dec_max = - fy_tyres(i) / M;
    ax_tyres_max = fx_tyres(i) / M;

    % Max Lateral acceleration
    ax = ax_tyres_max*cosd(linspace(0,180,ep))';

    % Max Longitudnal acceleration (Traction or Power Limited)
    ay_tyres = ay_tyres_acc_max * sqrt(1-(ax/ax_tyres_max).^2);
    ay_acc = (min(ay_tyres, ay_engine)) + ay_resistive;
    ay_dec = ay_tyres_dec_max * sqrt(1-(ax/ax_tyres_max).^2) + ay_resistive;
    
    % Generating GGV Map
    GGV(i,:,1) = [ay_acc',ay_dec(2:end)'] ;
    GGV(i,:,2) = [ax',flipud(ax(2:end))'] ;
    GGV(i,:,3) = vehicle_speed(i)*ones(1,2*ep-1) ;
end

disp('GGV map generated....')

w_last = warning('query','last');
warning('off', w_last.identifier);

%% ________________________________________________________________________
%% PLOTTING

% Engine Data
figure('Name',vehicle_name,'NumberTitle','off');
sgtitle('Powertrain & Traction Model');

subplot(2,2,1)
hold on;
title('Engine Torque/Power Curve');
yyaxis left;
plot(engine_RPM, engine_torque, LineWidth=1);
ylabel('Torque (Nm)')
yyaxis right;
plot(engine_RPM, engine_power * 746, LineWidth=1);
ylabel('Power (hp)');
xlabel('Engine RPM (rev/min)')
hold off;
grid on;

subplot(2,2,2);
title('Gear Torque Output');
plot(gear_speeds * 3.6, gear_torque, LineWidth=1);
ylabel('Torque (Nm)');
xlabel('Vehicle Speed (km/h)');
legend('1st Gear', '2nd Gear', '3rd Gear', '4th Gear', '5th Gear', '6th Gear', ...
    '7th Gear', '8th Gear', '9th Gear', '10th Gear');
grid on;

subplot(2,2,3);
hold on;
title('Gear Shift Points');
yyaxis left;
plot(vehicle_speed * 3.6, engine_RPM_mod, LineWidth=1);
ylabel('Engine RPM (rev/min');
yyaxis right;
plot(vehicle_speed * 3.6, gears, LineWidth=1);
ylim([0, max(gears) + 1]);
xlim([0, max(vehicle_speed) * 3.6]);
ylabel('nGear');
xlabel('Vehicle Speed (km/h)')
hold off;
grid on;

subplot(2,2,4);
hold on;
title('Vehicle Max. Traction Model');
plot(vehicle_speed * 3.6, fy_tyres, LineWidth=1);
plot(vehicle_speed * 3.6, tf_max, LineWidth=1.5, Color='k');
plot(vehicle_speed * 3.6, fy_aero, LineWidth=1);
plot(vehicle_speed * 3.6, gear_tf, LineStyle='-.', LineWidth=0.7, Color='k');
xlim([0, max(vehicle_speed)] * 3.6);
ylabel('Force (N)');
xlabel('Vehicle Speed (km/h)');
legend('Max. Tyre traction', 'Max. Engine traction', 'Drag', 'Engine Traction per Gear');
grid on;

% GGV Map
figure('Name',vehicle_name,'NumberTitle','off');
sgtitle('GGV Mapping');
surf(GGV(:,:,1)/g, GGV(:,:,2)/g, GGV(:,:,3) * 3.6, EdgeColor="none");
xlabel('Longitudnal Acceleration (g)');
ylabel('Lateral Acceleration (g)');
zlabel('Vehicle Speed (m/s)');
grid on;


%% ________________________________________________________________________
%% FUNCTIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[data] = read_specs(directory_name, excelFile, sheetName)
       
    % Using first sheet if sheet name is empty
    if nargin == 2 || isempty(sheetName)
        sheetName = 1;
    end

    startRow = 2;
    endRow = 10000;
    % Setup the Import Options
    opts = spreadsheetImportOptions("NumVariables", 2);
    % Specify sheet and range
    opts.Sheet = sheetName;
    opts.DataRange = "B" + startRow(1) + ":C" + endRow(1);
    % Specify column names and types
    opts.VariableNames = ["Spec", "Value"];
    opts.VariableTypes = ["string", "string"];
    % Setup rules for import
    opts.MissingRule = "omitrow";
    opts = setvaropts(opts, [1, 2], "TreatAsMissing", '');
    % Import the data
    data = readtable(append(directory_name, excelFile), opts, "UseExcel", false);
    for i = 2:length(startRow)
        opts.DataRange = "A" + startRow(i) + ":B" + endRow(i);
        tb = readtable(append(directory_name, excelFile), opts, "UseExcel", false);
        data = [data; tb];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data] = read_torque_curve(directory_name, excelFile, sheetName)
    startRow = 2;
    endRow = 10000;
    % Setup the Import Options
    opts = spreadsheetImportOptions("NumVariables", 2);
    % Specify sheet and range
    opts.Sheet = sheetName;
    opts.DataRange = "A" + startRow(1) + ":B" + endRow(1);
    % Specify column names and types
    opts.VariableNames = ["RPM", "Torque"];
    opts.VariableTypes = ["double", "double"];
    % Setup rules for import
    opts.MissingRule = "omitrow";
    opts = setvaropts(opts, [1, 2], "TreatAsMissing", '');
    % Import the data
    data = readtable(append(directory_name, excelFile), opts, "UseExcel", false);
    for i = 2:length(startRow)
        opts.DataRange = "A" + startRow(i) + ":B" + endRow(i);
        tb = readtable(append(directory_name, excelFile), opts, "UseExcel", false);
        data = [data; tb];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%