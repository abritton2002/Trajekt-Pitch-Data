% Clear the workspace and the command window
clear;
clc;

% Load the master dataset
master_file = 'Master_Data_Pitching_MLB.csv'; % Update this to the path of your CSV file
master_data = readtable(master_file);

% Load the active spin dataset
active_spin_file = 'Active-Spin.csv'; % Update this to the path of your CSV file
active_spin_data = readtable(active_spin_file);

% Load the pitch movement dataset
pitch_movement_file = 'pitch_movement.csv'; % Update this to the path of your CSV file
pitch_movement_data = readtable(pitch_movement_file);

% Convert active spin percentages to decimals and handle missing values
pitches = {'active_spin_fourseam', 'active_spin_sinker', 'active_spin_cutter', ...
           'active_spin_changeup', 'active_spin_splitter', 'active_spin_curve', ...
           'active_spin_slider', 'active_spin_sweeper'};
for pitch = pitches
    pitch_name = pitch{1};
    if iscell(active_spin_data.(pitch_name))
        % Convert cell array with percentage strings to numeric
        active_spin_data.(pitch_name) = str2double(strrep(active_spin_data.(pitch_name), '%', '')) / 100;
    else
        % Convert numeric percentages to decimals
        active_spin_data.(pitch_name) = active_spin_data.(pitch_name) / 100;
    end
end

% Extract the specified columns from the master data
columns_to_extract = {'player_name', 'pitch_type', 'spin_axis', 'release_spin_rate', 'vx0', 'vy0', 'vz0', 'p_throws', 'release_speed', 'release_pos_x', 'release_pos_z'};
extracted_data = master_data(:, columns_to_extract);

% Remove rows with any NaN values in the master data
extracted_data = rmmissing(extracted_data);

% Filter for specific pitch types
valid_pitch_types = {'FF', 'SI', 'FC', 'CH', 'FS', 'CU', 'SL', 'ST'};
extracted_data = extracted_data(ismember(extracted_data.pitch_type, valid_pitch_types), :);

% Get unique combinations of player_name and pitch_type
unique_combinations = unique(extracted_data(:, {'player_name', 'pitch_type'}));

% Initialize an empty table for storing results
results = table();

% Iterate through each unique combination of player_name and pitch_type
for i = 1:height(unique_combinations)
    player_name = unique_combinations.player_name{i};
    pitch_type = unique_combinations.pitch_type{i};

    % Filter data for the current player_name and pitch_type
    filtered_data = extracted_data(strcmp(extracted_data.player_name, player_name) & ...
                                   strcmp(extracted_data.pitch_type, pitch_type), :);

    % Extract the pitcher's throwing hand (left or right)
    p_throws = filtered_data.p_throws{1}; % Extract the first cell as a string

    % Calculate the average values
    avg_spin_axis = mean(filtered_data.spin_axis);
    avg_release_spin_rate = mean(filtered_data.release_spin_rate);
    avg_vx0 = mean(filtered_data.vx0);
    avg_vy0 = mean(filtered_data.vy0);
    avg_vz0 = mean(filtered_data.vz0);
    avg_release_speed = mean(filtered_data.release_speed);
    avg_release_pos_x = mean(filtered_data.release_pos_x);
    avg_release_pos_z = mean(filtered_data.release_pos_z);

    % Split the player name into first and last names
    name_parts = split(player_name, ',');
    last_name = strtrim(name_parts{1});
    first_name = strtrim(name_parts{2});
    pitcher_full_name = [first_name, ' ', last_name];

    % Find the spin efficiency for the current player and pitch type
    spin_efficiency_row = active_spin_data(strcmp(active_spin_data.('last_name_First_name'), [last_name, ', ', first_name]) & ...
                                           strcmp(active_spin_data.pitch_hand, p_throws), :);

    if isempty(spin_efficiency_row)
        fprintf('No spin efficiency data found for %s (%s).\n', player_name, pitch_type);
        continue;
    end

    % Extract the spin efficiency value based on the pitch type
    switch pitch_type
        case 'FF'
            spin_efficiency = spin_efficiency_row.active_spin_fourseam;
        case 'SI'
            spin_efficiency = spin_efficiency_row.active_spin_sinker;
        case 'FC'
            spin_efficiency = spin_efficiency_row.active_spin_cutter;
        case 'CH'
            spin_efficiency = spin_efficiency_row.active_spin_changeup;
        case 'FS'
            spin_efficiency = spin_efficiency_row.active_spin_splitter;
        case 'CU'
            spin_efficiency = spin_efficiency_row.active_spin_curve;
        case 'SL'
            spin_efficiency = spin_efficiency_row.active_spin_slider;
        case 'ST'
            spin_efficiency = spin_efficiency_row.active_spin_sweeper;
        otherwise
            spin_efficiency = NaN;
    end

    % Ensure spin_efficiency is a scalar
    spin_efficiency = spin_efficiency(1);

    % If spin efficiency is not found or NaN, skip this entry
    if isnan(spin_efficiency)
        fprintf('No spin efficiency data found for %s (%s).\n', player_name, pitch_type);
        continue;
    end

    % Calculate the average pitch movement values for the pitcher and pitch type
    pitcher_movement_data = pitch_movement_data(strcmp(pitch_movement_data.('last_name_First_name'), [last_name, ', ', first_name]) & ...
                                                strcmp(pitch_movement_data.pitch_type, pitch_type), :);
    if ~isempty(pitcher_movement_data)
        avg_pitcher_break_x = mean(pitcher_movement_data.pitcher_break_x);
        avg_pitcher_break_z = mean(pitcher_movement_data.pitcher_break_z);

        % Adjust the break_x value if the pitcher is a righty
        if strcmp(p_throws, 'R')
            avg_pitcher_break_x = -avg_pitcher_break_x;
        end
    else
        avg_pitcher_break_x = NaN;
        avg_pitcher_break_z = NaN;
        fprintf('No pitch movement data found for %s (%s).\n', player_name, pitch_type);
    end

    % Total spin rate (in rpm)
    omega = avg_release_spin_rate;

    % Convert spin_axis to radians
    phi = deg2rad(avg_spin_axis);

    % Convert velocities from feet/s to meters/s
    feet_to_meters = 0.3048;
    vx = avg_vx0 * feet_to_meters;
    vy = avg_vy0 * feet_to_meters;
    vz = avg_vz0 * feet_to_meters;

    % Total velocity
    v = sqrt(vx^2 + vy^2 + vz^2);

    % Calculate the horizontal (ϕ) and vertical (θ) release angles
    phi_release = atan2(vx, vy);
    theta_release = asin(vx / v);

    % Velocity unit vectors using release angles
    v_hat_x = cos(theta_release) * sin(phi_release);
    v_hat_y = cos(theta_release) * cos(phi_release);
    v_hat_z = sin(theta_release);

    % Calculate cos(theta_S)
    cos_theta_S = sqrt(1 - spin_efficiency.^2);

    % Calculate parameters A, B, C
    A = v_hat_x * cos(phi) + v_hat_z * sin(phi);
    B = v_hat_y;
    C = cos_theta_S;

    % Calculate R and X
    R = sqrt(A^2 + B^2);
    X = atan2(B, A);

    % Calculate the polar angle Theta
    Theta = asin(C / R) - X;

    % Calculate the spin components
    omega_x = omega * sin(Theta) * cos(phi);
    if strcmp(p_throws, 'L')
        omega_y = -1 * omega * cos(Theta);
    else
        omega_y = omega * cos(Theta);
    end
    omega_z = omega * sin(Theta) * sin(phi);

    % Store the results in the table
    new_row = table({pitcher_full_name}, {pitch_type}, ...
                    avg_release_pos_x, 56.5, avg_release_pos_z, avg_release_speed, ...
                    omega_x, omega_y, omega_z, avg_pitcher_break_x, avg_pitcher_break_z, ...
                    'VariableNames', {'Player_Name', 'PitchType', ...
                    'Release_Pos_X', 'ReleaseY', 'Release_Pos_Z', 'Release_Speed', ...
                    'Spin_X', 'Spin_Y', 'Spin_Z', 'Break_X', 'Break_Z'});
    results = [results; new_row];
end

% Define the Excel file name
excel_file = 'pitch_data_master_v4_MLB.xlsx';

% Write the results to the Excel file
writetable(results, excel_file);

fprintf('Data successfully saved to %s\n', excel_file);
