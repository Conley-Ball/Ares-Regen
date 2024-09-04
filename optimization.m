clear all; close all; clc;

% start and end values
w_ch_min_start = 0.001; % m
w_rib_start = 0.001; % m
t_ins_start = 0.001; % m
h_ch_th_start = 0.001; % m
h_ch_c_start = 0.001; % m
h_ch_e_start = 0.001; % m

h_ch_th_end = 0.001; % m
h_ch_c_end = 0.003; % m
h_ch_e_end = 0.003; % m
w_ch_min_end = 0.001; % m
w_rib_end = 0.003; % m
t_ins_end = 0.0035; % m

% number of values
num_vals_w_ch_min = 1;
num_vals_w_rib = 5;
num_vals_t_ins = 5;
num_vals_h_ch_th = 1;
num_vals_h_ch_c = 5;
num_vals_h_ch_e = 5;

% CEA Inputs
Pc = 313.7; % psia
Pe = 13.7; % psia
O_F = .9;
T_inlet = 300; % K
res = 0.00005/1; % thermal resistance coating
Thrust = 2000; % lbf
C_star_eff = 0.94;
C_F_eff = 0.99;
L_star = 30; % in
t_out = 0.001 / 0.0254; % in
ratio = 0.75;

% CEA Function
[AR, C_star, C_F, gamma, MW_g, Cp_g, mu_g, k_g, T_thr, Pr_g] = runCEA(Pc, Pe, O_F, ratio);
T_thr = T_thr * C_star_eff^2;

% generate values
w_ch_min = linspace(w_ch_min_start, w_ch_min_end, num_vals_w_ch_min);
w_rib = linspace(w_rib_start, w_rib_end, num_vals_w_rib);
t_ins = linspace(t_ins_start, t_ins_end, num_vals_t_ins);
h_ch_th = linspace(h_ch_th_start, h_ch_th_end, num_vals_h_ch_th);
h_ch_c = linspace(h_ch_c_start, h_ch_c_end, num_vals_h_ch_c);
h_ch_e = linspace(h_ch_e_start, h_ch_e_end, num_vals_h_ch_e);

% convert to inches
h_ch_th = h_ch_th / 0.0254; % in
h_ch_c = h_ch_c / 0.0254; % in
h_ch_e = h_ch_e / 0.0254; % in
w_ch_min = w_ch_min / 0.0254; % in
t_ins = t_ins / 0.0254; % in
w_rib = w_rib / 0.0254; % in

% Determine number of cases
num_cases = numel(w_ch_min) * numel(w_rib) * numel(t_ins) * numel(h_ch_th) * numel(h_ch_c) * numel(h_ch_e);

% Initialize variables to track FOS
max_fos = -inf; % Initialize to negative infinity to find the maximum
best_case = -1;

% initialize results
results = cell(num_cases, 11); 

% counter
case_index = 1;

% loop through all parameter combinations
for i = 1:numel(w_ch_min)
    for j = 1:numel(w_rib)
        for k = 1:numel(t_ins)
            for l = 1:numel(h_ch_th)
                for m = 1:numel(h_ch_c)
                    for n = 1:numel(h_ch_e)
                        % get current parameters
                        current_w_ch_min = w_ch_min(i);
                        current_w_rib = w_rib(j);
                        current_t_ins = t_ins(k);
                        current_h_ch_th = h_ch_th(l);
                        current_h_ch_c = h_ch_c(m);
                        current_h_ch_e = h_ch_e(n);

                        % call the main function
                        [v_m_stress, Yield_iw, T_chg, P_c] = main_function(current_w_ch_min, current_w_rib, current_t_ins, current_h_ch_th, current_h_ch_c, current_h_ch_e, Pc, Pe, O_F, T_inlet, res, Thrust, C_star_eff, C_F_eff, L_star, t_out, ratio, AR, C_star, C_F, gamma, MW_g, Cp_g, mu_g, k_g, T_thr, Pr_g);

                        % calculate FOS for each position
                        fos_array = Yield_iw ./ v_m_stress;

                        % find the minimum FOS for this parameter combination
                        min_fos = min(fos_array);

                        % store results for the minimum FOS
                        results{case_index, 1} = current_w_ch_min * 0.0254; % meters
                        results{case_index, 2} = current_w_rib * 0.0254; % meters
                        results{case_index, 3} = current_t_ins * 0.0254; % meters
                        results{case_index, 4} = current_h_ch_th * 0.0254; % meters
                        results{case_index, 5} = current_h_ch_c * 0.0254; % meters
                        results{case_index, 6} = current_h_ch_e * 0.0254; % meters
                        results{case_index, 7} = min(v_m_stress) / 6.895e+6; % ksi
                        results{case_index, 8} = min(Yield_iw) / 6.895e+6; % ksi
                        results{case_index, 9} = max(T_chg); % K
                        results{case_index, 10} = min_fos; % Minimum FOS
                        dP = P_c(1) - P_c(end);
                        results{case_index, 11} = dP;

                        % check for highest minimum FOS and dp<120ksi
                        if min_fos > max_fos && dP < 120
                            max_fos = min_fos;
                            best_case = case_index;
                        end

                        % increase index
                        case_index = case_index + 1;
                    end
                end
            end
        end
    end
end

% generate table
results_table = cell2table(results, 'VariableNames', {'w_ch_min_m', 'w_rib_m', 't_ins_m', 'h_ch_th_m', 'h_ch_c_m', 'h_ch_e_m', 'v_m_stress_ksi', 'Yield_iw_ksi', 'T_chg_K', 'FOS', 'dP'});

% print the case with the highest minimum FOS
if best_case ~= -1
    best_case_table = results_table(best_case, :);
    fprintf('Case with highest minimum FOS:\n');
    disp(best_case_table);
else
    fprintf('No cases found.\n');
end
