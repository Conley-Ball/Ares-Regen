clear all; close all; clc;

% start values (don't decrease below min values)
w_ch_min_start = 0.000635; % 0.000635 min
w_rib_start = 0.000381; % 0.000381 min
t_ins_start = 0.00042; % 0.00042 min
h_ch_start = 0.000635; % 0.000635 min

% end values
w_ch_min_end = 0.000636;
w_rib_end = 0.000382;
t_ins_end = 0.00043;
h_ch_end = 0.000635;

% number of values for each parameter
num_vals_w_ch_min = 2;
num_vals_w_rib = 2;
num_vals_t_ins = 2;
num_vals_h_ch = 1;

%CEA Inputs
Pc = 413; % psia
Pe = 13.7; % psia
O_F = 1.1;
T_inlet = 300; % K

res = 0; % thermal resistance coating
Thrust = 1600; % lbf
%Pe = 13; % psi 
C_star_eff = 0.94;
C_F_eff = 0.9;
L_star = 28; % in

t_out = 0.001/0.0254; % in

ratio = 0.75;

%% CEA Function

[AR, C_star, C_F, gamma, MW_g, Cp_g, mu_g, k_g, T_thr, Pr_g] = runCEA(Pc, Pe, O_F, ratio);
T_thr = T_thr*C_star_eff;

%%
% generate values
w_ch_min = linspace(w_ch_min_start, w_ch_min_end, num_vals_w_ch_min);
w_rib = linspace(w_rib_start, w_rib_end, num_vals_w_rib);
t_ins = linspace(t_ins_start, t_ins_end, num_vals_t_ins);
h_ch = linspace(h_ch_start, h_ch_end, num_vals_h_ch);

%convert to inches
h_ch = h_ch/0.0254; % in
w_ch_min = w_ch_min/0.0254; % in
t_ins = t_ins/0.0254; % in
w_rib = w_rib/0.0254; % in


% determine number of cases
num_cases = numel(w_ch_min) * numel(w_rib) * numel(t_ins) * numel(h_ch);
% initialize results
results = cell(num_cases, 8);  

% initialize variables to track FOS
max_fos = 0;  
best_case = -1;   
% counter
case_index = 1;
% Loop through all parameter combinations
for i = 1:numel(w_ch_min)
    for j = 1:numel(w_rib)
        for k = 1:numel(t_ins)
            for l = 1:numel(h_ch)
                % get current parameters
                current_w_ch_min = w_ch_min(i);
                current_w_rib = w_rib(j);
                current_t_ins = t_ins(k);
                current_h_ch = h_ch(l);
                % call the main function
                [v_m_stress, Yield_iw, T_chg] = main_function(current_w_ch_min, current_w_rib, current_t_ins, current_h_ch, Pc,Pe,O_F, T_inlet,res,Thrust,C_star_eff,C_F_eff,L_star,t_out,ratio,AR, C_star, C_F, gamma, MW_g, Cp_g, mu_g, k_g, T_thr, Pr_g);

                
                % determine the relevant values
                max_v_m_stress = max(v_m_stress);
                min_Yield_iw = min(Yield_iw);
                max_T_chg = max(T_chg);
                % determine fos
                fos = min_Yield_iw / max_v_m_stress;
                % store results
                results{case_index, 1} = current_w_ch_min;
                results{case_index, 2} = current_w_rib;
                results{case_index, 3} = current_t_ins;
                results{case_index, 4} = current_h_ch;
                results{case_index, 5} = max_v_m_stress/6.895e+6;
                results{case_index, 6} = min_Yield_iw/6.895e+6;
                results{case_index, 7} = max_T_chg;
                results{case_index, 8} = fos;
                
                % check if this is the highest FOS so far
                if fos > max_fos
                    max_fos = fos;
                    best_case = case_index;
                end
                
                % increase index
                case_index = case_index + 1;
            end
        end
    end
end

% generate table
resutls_table = cell2table(results, 'VariableNames', {'w_ch_min (in)', 'w_rib (in)', 't_ins (in)', 'h_ch (in)', 'v_m_stress (ksi)', 'Yield_iw (ksi)', 'T_chg (k)', 'FOS'});

% print the case with the highest FOS
if best_case ~= -1
    best_case_table = resutls_table(best_case, :);
    fprintf('Case with highest FOS:\n');
    disp(best_case_table);
else
    fprintf('No cases found.\n');
end