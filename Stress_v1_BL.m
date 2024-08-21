function [stress_total,stress_T_ow, stress_T_iw, stress_T_s, stress_T_c,stress_P_b, stress_P_c, stress_P_s, stress_P_t, stress_P_hoop, stress_P_a] = stress(T_iw,T_ow,E_iw,E_ow,E_avg,alpha_iw,alpha_ow,alpha_avg,nu_iw,nu_ow,nu_avg,t_iw,t_ow,h,w,d,r_ow,r_iw,deltaP,P_chamber,A_wall,A_cc,numel)
%{ =====STRESS===== 
% Total Stress calculation at each cross section

% Sources: 
% (1) DEVELOPMENT OF ADVANCED FABRICATION TECHNIQUES FOR REGENERATIVELY COOLED THRUST CHAMBERS
%}

% Define arrays
[stress_T_ow, stress_T_iw, stress_T_s, stress_T_c] = deal(zeros(numel,1));
[stress_P_b, stress_P_c, stress_P_s, stress_P_t, stress_P_hoop, stress_P_a] = deal(zeros(numel,1));

% 
for i=1:numel

% THERMAL STRESS [Pa]
% Circumferential stress outer wall  [(1)B3,p120] 
stress_T_ow(i) = ((alpha_ow(i)*E_ow(i)*(T_iw(i)-T_ow(i)))/(1-nu_ow(i)))*(t_iw(i)/(t_iw(i)+t_ow(i)));
% Circumferential stress inner wall [(1)B4,p121] 
stress_T_iw(i) = ((alpha_iw(i)*E_iw(i)*(T_iw(i)-T_ow(i)))/(1-nu_iw(i)))*(t_ow(i)/(t_iw(i)+t_ow(i)));
% Shear Stress max [(1)B4,p121]
stress_T_s(i) = ( (alpha_avg(i)*E_avg(i)*(T_iw(i)-T_ow(i)))/(5*w(i)*(1-nu_avg(i))) ) * ((t_ow(i)*t_iw(i))/(t_iw(i)+t_ow(i))^2) * (h(i)+w(i));
% Compression Stress max [(1)B4,p121]
stress_T_c(i) = ( (alpha_avg(i)*E_avg(i)*(T_iw(i)-T_ow(i))) / ((w(i)*r_ow(i)/(h(i)+w(i)))*(1-nu_avg(i))) ) * ((t_ow(i)*t_iw(i))/(t_iw(i)+t_ow(i)));

% PRESSURE STRESS: INNER WALL [Pa]
% Maximum bending stress  [(1)A3,p111]
stress_P_b(i) = (deltaP/2)*(h(i)/t_iw(i))^2;
% Compressive loading, ribs
stress_P_c(i) = (deltaP/w(i))*(h(i)+w(i));
% Tensile stress, ribs
stress_P_t(i) = (deltaP/w(i))*h(i);
% Shear stress, ribs
stress_P_s(i) = (deltaP/2)*(h(i)/t_iw(i));
% circumferential stress in the inner and outer walls
stress_P_hoop(i) = deltaP*r_ow(i)/(t_iw(i)+t_ow(i));
% Hoop 
stress_P_hoop(i) = deltaP*r_iw(i)/(t_iw(i)+t_ow(i));
% Axial 
stress_P_a(i) = P_chamber.*A_cc./A_wall(i);

end


% TOTAL STRESS
% these stresses can be simply added together to determine the total stress on an object because of the superposition principle.
stress_total = stress_T_ow + stress_T_iw + stress_T_s + stress_T_c + stress_P_b + stress_P_c + stress_P_s + stress_P_t + stress_P_hoop + stress_P_a;

% Von Mises maximum distortion energy theory
% failure by yielding occurs when, at any point in the body, the distortion energy per unit volume 
% in a state of combined stress becomes equal to that associated with yielding

end

