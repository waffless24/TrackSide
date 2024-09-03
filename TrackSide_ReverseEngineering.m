%% ________________________________________________________________________
%% TRACK SIDE VEHICLE MODELLING SCRIPT
%
% Vehicle Model generation for use in lap-time simulation
%
% Liscensing: GPL V3 Open Source License.
% Information:
%   Nirat Pai
%   BSc Mechanical Engineering, Purdue University
%   LinkedIn: linkedin.com/in/nirat-pai-770661288
%   email: niratpai24@gmail.com
%   GitHub: https://github.com/waffless24

%% ________________________________________________________________________
%% CLEARING WORKSPACE
clear
clc

%% ________________________________________________________________________
%% GENERIC PARAMETERS

% TYRES
tyre_rad = 360/1000;

% TRANSMISSION
r_primary = 1;
r_drive = 4.9;
e_driveline = 0.9061;
nog = 8;
r_gear = zeros(1, nog).';



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
    vehicle_speed = telemetry{:,"Speed"} * (5/18);
    engine_RPM = telemetry{:, "RPM"};
    ngear = telemetry{:, "nGear"};

%% ________________________________________________________________________
%% GEAR RATIOS
    
    r_final = ((engine_RPM * (pi/30)) ./ vehicle_speed) * tyre_rad;
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
    

    

end

%% ________________________________________________________________________
%% FINAL RESULTS

for i = 1:nog
    r_gear(i) = mean(r_gear(i,:), "omitmissing");
end

r_gear(:,2:end) = [];

