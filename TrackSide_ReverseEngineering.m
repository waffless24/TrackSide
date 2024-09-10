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
%% GENERIC PARAMETERS (ASSUMED)

g = 9.81;
fuel_M = 6.132; %[kg]

% TYRES
directory = dir('Vehicle Parts/Tyres & Steering/*.xlsx');
c_soft = directory(find(contains({directory.name},'F1 Slicks Soft'))).name;
c_soft = read_specs('Vehicle Parts/Tyres & Steering/', c_soft);

c_medium = directory(find(contains({directory.name},'F1 Slicks Medium'))).name;
c_medium = read_specs('Vehicle Parts/Tyres & Steering/', c_medium);

c_hard = directory(find(contains({directory.name},'F1 Slicks Hard'))).name;
c_hard = read_specs('Vehicle Parts/Tyres & Steering/', c_hard);

tyre_rad = (double(c_soft{3, "Value"}) * 0.0254) / 2;

% Coefficient of Rolling
mu_r = [double(c_soft{4, "Value"}), double(c_medium{4, "Value"}), double(c_hard{4, "Value"})];

% Longitudnal Coefficients
mu_y = [double(c_soft{5, "Value"}), double(c_medium{5, "Value"}), double(c_hard{5, "Value"})];
loadR_y = [double(c_soft{6, "Value"}), double(c_medium{6, "Value"}), double(c_hard{6, "Value"})];
mu_y_sens = [double(c_soft{7, "Value"}), double(c_medium{7, "Value"}), double(c_hard{7, "Value"})];

% Lateral Coefficients
mu_x = [double(c_soft{8, "Value"}), double(c_medium{8, "Value"}), double(c_hard{8, "Value"})];
loadR_x = [double(c_soft{9, "Value"}), double(c_medium{9, "Value"}), double(c_hard{9, "Value"})];
mu_x_sens = [double(c_soft{10, "Value"}), double(c_medium{10, "Value"}), double(c_hard{10, "Value"})];

% Tyre stiffness (Torque Resistance)
stiff_F = [double(c_soft{11, "Value"}), double(c_medium{11, "Value"}), double(c_hard{11, "Value"})];
stiff_R = [double(c_soft{12, "Value"}), double(c_medium{12, "Value"}), double(c_hard{12, "Value"})];

% TRANSMISSION
r_primary = 1;
r_drive = 4.9;
e_driveline = 0.9061;
nog = 8;
r_gear = zeros(1, nog).';

% CHASSIS
M = 765;
M_dist = 0.46;
aero_dist = 0.55;

%% ________________________________________________________________________
%% READING FILES

% Set the desired vehicle
vehicle_name = 'RB20 Package 1';

% Obtaining files from folder
directory_name = append('Reverse Engineering Data/', vehicle_name, '/');
directory = dir(append('Reverse Engineering Data/', vehicle_name, '/*.csv'));

% Starting loop to read through each file
for i = 1:length(directory)
    filename = directory(i).name;
    telemetry = readtable(append(directory_name, filename));
    c_index = 1;
    vehicle_speed = telemetry{:,"Speed"} * (5/18);
    engine_RPM = telemetry{:, "RPM"};
    ngear = telemetry{:, "nGear"};

%% ________________________________________________________________________
%% RACING LINE PARAMETERS
    
    % Resetting previous entries
    X = [];
    Y = [];
    Z = [];
    radius = [];
    incl = [];
    bank = [];
    
    X = telemetry{:, "X"} * 0.1;
    Y = telemetry{:, "Y"} * 0.1;
    Zmin = min(telemetry{:, "Z"});
    Z = (telemetry{:, "Z"} - Zmin) * 0.1;

    X = X - X(1);
    Y = Y - Y(1);
    Z = Z - Z(1);

    for j = 1:height(X)      
        if j+2 > height(X)
            if j+1 > height(X)
                A = [X(j), Y(j), Z(j)];
                B = [X(1), Y(1), Z(1)];
                C = [X(2), Y(2), Z(2)];
            else
                A = [X(j), Y(j), Z(j)];
                B = [X(j+1), Y(j+1), Z(j+1)];
                C = [X(1), Y(1), Z(1)];
            end
        else
            A = [X(j), Y(j), Z(j)];
            B = [X(j+1), Y(j+1), Z(j+1)];
            C = [X(j+2), Y(j+2), Z(j+2)];
        end
    
    [radius(j,1), incl(j,1), bank(j,1)] = get_radius(A,B,C);
    end

%% ________________________________________________________________________
%% GEAR RATIOS
     
    r_final = ((engine_RPM * (pi/30)) ./ vehicle_speed) * tyre_rad(c_index);
    r_gears  = r_final / (r_primary * r_drive);

    r_gears = [r_gears, ngear];

    for j = 1:nog
         if ismember(j, r_gears(:,2))
            gear = r_gears(ismember(r_gears(:,2),j),:);

            % Filtering out outliers
            gear = rmoutliers(gear,"median");

            r_gear(j,i) = mean(gear(:,1));
         else
            r_gear(j,i) = NaN;
            continue
         end
    end

%% ________________________________________________________________________
%% AERO


end
%% ________________________________________________________________________
%% FINAL RESULTS

% TRANSMISSION
for i = 1:nog
    r_gear(i) = mean(r_gear(i,:), "omitmissing");
end

r_gear(:,2:end) = [];


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

function[radius, incl, bank] = get_radius(A, B, C)
       
    u1 = B - A;
    w1 = cross(C-A, u1);
    u = u1/norm(u1);
    w = w1/norm(w1);
    v = cross(w,u);

    b = [dot((B-A), u), 0];
    c = [dot((C-A), u), dot((C-A), v)];

    h1 = (c(1)-b(1)/2)^2 + c(2)^2 -(b(1)/2)^2;
    h = h1/(2 * c(2));

    centre(1,:) = A + (b(1)/2)*u + h*v;
    radius = norm(A - centre(1,:));
    
    incl = A(3)/radius;
    bank = 90 - acosd(1/radius);
    radius = radius * cosd(bank);

    if radius >= 500
        radius = NaN;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
