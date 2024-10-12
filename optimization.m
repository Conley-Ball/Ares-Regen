clear all; close all; clc;
% Manually adjust the number of values for each parameter cuz im tired and
% havent implemented it yet
%% Geometry Inputs
% Define channel heights (start and end) and convert from inches to meters
current_h_ch = [0 0 0 0];
% Define insulation thicknesses (start and end) and convert from inches to meters
t_ins_start = [0.002 0.001 0.001 0.001];% m
t_ins_end = [0.002 0.002 0.002 0.002];% m
t_ins_start = t_ins_start / 0.0254; % in
t_ins_end = t_ins_end / 0.0254; % in

% Generate ranges for insulation thicknesses
t_ins_1 = linspace(t_ins_start(1), t_ins_end(1), 1);% in
t_ins_2 = linspace(t_ins_start(2), t_ins_end(2), 2);% in
t_ins_3 = linspace(t_ins_start(3), t_ins_end(3), 2);% in
t_ins_4 = linspace(t_ins_start(4), t_ins_end(4), 2);% in

% Define minimum channel width and convert from inches to meters
w_ch_min_start = 0.001; % m
w_ch_min_end = 0.001; % m
w_ch_min_start = w_ch_min_start / 0.0254; % in
w_ch_min_end = w_ch_min_end / 0.0254; % in
w_ch_min = linspace(w_ch_min_start, w_ch_min_end, 1); % in

% Define rib width and convert from inches to meters
w_rib_start = 0.001; % m
w_rib_end = 0.001; % m
w_rib_start = w_rib_start / 0.0254; % in
w_rib_end = w_rib_end / 0.0254; % in
w_rib = linspace(w_rib_start, w_rib_end, 1);% in

% Define output thickness and convert from inches to meters
t_out_start = [0.001 0.001 0.001 0.001];% m
t_out_end = [0.002 0.002 0.002 0.002];% m
t_out_start = t_out_start / 0.0254; % in
t_out_end = t_out_end / 0.0254; % in

% Generate ranges for outer wall thicknesses
t_out_1 = linspace(t_out_start(1), t_out_end(1), 2);% in
t_out_2 = linspace(t_out_start(2), t_out_end(2), 2);% in
t_out_3 = linspace(t_out_start(3), t_out_end(3), 2);% in
t_out_4 = linspace(t_out_start(4), t_out_end(4), 2);% in

%%
% CEA Inputs
Pc = 313.7; % psia
Pe = 10.2; % psia
O_F = .9;
T_inlet = 300; % K
res = 0.00005 / 1; % thermal resistance coating
Thrust = 2000; % lbf
C_star_eff = 0.94;
C_F_eff = 0.99;
L_star = 30; % in
ratio = 0.75;

fillet = 0.0005; % m radius

% Call CEA function
[AR, C_star, C_F, gamma, MW_g, Cp_g, mu_g, k_g, T_thr, Pr_g] = runCEA(Pc, Pe, O_F, ratio);
T_thr = T_thr * C_star_eff^2;

% Initialize results
num_cases = numel(w_ch_min) * numel(w_rib) * numel(t_ins_1) * numel(t_ins_2) * numel(t_ins_3) * numel(t_ins_4) * numel(t_out_1) * numel(t_out_2) * numel(t_out_2) * numel(t_out_3);
results = cell(num_cases, 13);
max_fos = -inf;
best_case = -1;

% Initialize progress
case_index = 1;
h_waitbar = waitbar(0, 'Progress: 0% Complete');

% Loop through all combinations
for i = 1:numel(w_ch_min)
    for j = 1:numel(w_rib)
        for k = 1:numel(t_ins_1)
            for p = 1:numel(t_ins_2)
            for q = 1:numel(t_ins_3)
            for r = 1:numel(t_ins_4)
            for l = 1:numel(t_out_1)
                for m = 1:numel(t_out_2)
                    for n = 1:numel(t_out_3)
                        for o = 1:numel(t_out_4)
                            % Get current parameters
                            current_w_ch_min = w_ch_min(i);
                            current_w_rib = w_rib(j);
                            current_t_ins_1 = t_ins_1(k);
                            current_t_ins_2 = t_ins_2(p);
                            current_t_ins_3 = t_ins_3(q);
                            current_t_ins_4 = t_ins_4(r);
                            current_t_out_1 = t_out_1(l);
                            current_t_out_2 = t_out_2(m);
                            current_t_out_3 = t_out_3(n);
                            current_t_out_4 = t_out_4(o);
                            
                            current_t_ins = [current_t_ins_1,current_t_ins_2,current_t_ins_3,current_t_ins_4];
                            current_t_out = [current_t_out_1,current_t_out_2,current_t_out_3,current_t_out_4];

                            if max_fos == -inf %for first case/none found yet, base fos value
                                fos = 1.1;
                            else
                                fos = max_fos+0.01; %optimization, no point in searching for fos under max fos
                            end
                            valid = true; %false if min feature size is reached in heatTransfer2D
                            while valid
                            % Call main
                            [T_chg, P_c, valid] = main_function(current_w_ch_min, current_w_rib, current_t_ins, current_h_ch, Pc, Pe, O_F, T_inlet, res, Thrust, C_star_eff, C_F_eff, L_star, current_t_out, ratio, AR, C_star, C_F, gamma, MW_g, Cp_g, mu_g, k_g, T_thr, Pr_g, fos, valid, fillet);
                            % Store results
                            results{case_index, 1} = current_w_ch_min * 0.0254; % m
                            results{case_index, 2} = current_w_rib * 0.0254; % m
                            results{case_index, 3} = current_t_ins_1 * 0.0254; % m
                            results{case_index, 4} = current_t_ins_2 * 0.0254; % m
                            results{case_index, 5} = current_t_ins_3 * 0.0254; % m
                            results{case_index, 6} = current_t_ins_4 * 0.0254; % m

                            results{case_index, 7} = current_t_out_1 * 0.0254; % m
                            results{case_index, 8} = current_t_out_2 * 0.0254; % m
                            results{case_index, 9} = current_t_out_3 * 0.0254; % m
                            results{case_index, 10} = current_t_out_4 * 0.0254; % m
                            results{case_index, 11} = max(T_chg); % K
                            results{case_index, 12} = fos; % Min FOS
                            dP = P_c(1) - P_c(end);
                            results{case_index, 13} = dP;

                            %test if this geometry can have higher fos
                            fos = fos+0.01;
                            end
                            if fos == 1.11 %if the geometry did not pass the initial fos
                                fos = -inf;
                            else
                                fos = fos-0.02; %if the geometry passed the initial fos, must subtract due to order of while loop
                            end
                            % Check for highest minimum FOS
                            if fos > max_fos 
                                max_fos = fos;
                                best_case = case_index;
                            end

                            % Increase index
                            case_index = case_index + 1;
                            % Track progress
                            progress = (i-1) * numel(w_rib) * numel(t_ins_1) * numel(t_ins_2) * numel(t_ins_3) * numel(t_ins_4) * numel(t_out_1) * numel(t_out_2) * numel(t_out_3) * numel(t_out_4) + ...
                                                   (j-1) * numel(t_ins_1) * numel(t_ins_2) * numel(t_ins_3) * numel(t_ins_4) * numel(t_out_1) * numel(t_out_2) * numel(t_out_3) * numel(t_out_4) + ...
                                                   (k-1) * numel(t_ins_2) * numel(t_ins_3) * numel(t_ins_4) * numel(t_out_1) * numel(t_out_2) * numel(t_out_3) * numel(t_out_4) + ...
                                                   (p-1) * numel(t_ins_3) * numel(t_ins_4) * numel(t_out_1) * numel(t_out_2) * numel(t_out_3) * numel(t_out_4) + ...
                                                   (q-1) * numel(t_ins_4) * numel(t_out_1) * numel(t_out_2) * numel(t_out_3) * numel(t_out_4) + ...
                                                   (r-1) * numel(t_out_1) * numel(t_out_2) * numel(t_out_3) * numel(t_out_4) + ...
                                                   (l-1) * numel(t_out_2) * numel(t_out_3) * numel(t_out_4) + ...
                                                   (m-1) * numel(t_out_3) * numel(t_out_4) + ...
                                                   (n-1) * numel(t_out_4) + ...
                                                   o;
                                        waitbar(progress / num_cases, h_waitbar, sprintf('Progress: %.2f%% Complete', (progress / num_cases) * 100));
                        end
                    end
                end
            end
            end
            end
            end
        end
    end
end
close(h_waitbar);

% Generate table
results_table = cell2table(results, 'VariableNames', {'w_ch_min_m', 'w_rib_m', 't_ins_1', 't_ins_2', 't_ins_3', 't_ins_4', 't_out_1', 't_out_2', 't_out_3', 't_out_4','T_chg_K', 'FOS', 'dP'});

% Print the case with the highest minimum FOS
if best_case ~= -1
    best_case_table = results_table(best_case, :);
    fprintf('Case with highest minimum FOS:\n');
    disp(best_case_table);
else
    fprintf('No cases')
end