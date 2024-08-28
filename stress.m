function [stressTow, stressTiw, stressTs, stressTc, stressPb, stressPs, stressPt, stressP_hoop, stressPa, vonMises] = stress(T_chg,T_ro,E_iw,E_ow,alpha_iw,alpha_ow,nu_iw,nu_ow,t_ins,t_out,w_ch,w_rib,h_ch,D_c,P_c,P,num_ch,A_t)
%{ =====STRESS===== 
% Total Stress calculation at each cross section

% Sources: 
% (1) DEVELOPMENT OF ADVANCED FABRICATION TECHNIQUES FOR REGENERATIVELY COOLED THRUST CHAMBERS
%}

% Define arrays
% [stress_T_ow, stress_T_iw, stress_T_s, stress_T_c] = deal(zeros(1,numel));
% [stress_P_b, stress_P_c, stress_P_s, stress_P_t, stress_P_hoop, stress_P_a] = deal(zeros(1,numel));

% Average inner and outer wall
E_avg = (E_iw+E_ow)/2;
alpha_avg = (alpha_iw+alpha_ow)/2;
nu_avg = (nu_iw+nu_ow)/2;

deltaP = P_c-P;  % Pressure difference [Pa] channel(higher) to chamber (30  psi)
P_chamber = P;

% Engine and Channel dimensions 
r_iw = D_c/2; 
r_ow = r_iw+t_ins+h_ch+t_out;
r_bar = (r_iw+r_ow)/2;
A_channels = num_ch*w_ch*h_ch;
A_wall = pi*(r_ow.^2 - r_iw.^2) - A_channels;
A_chamber = pi*r_iw.^2;
A_cc = A_chamber-A_t;


% THERMAL STRESS [Pa]
% Circumferential stress outer wall  [(1)B3,p120] ()
stressTow = ((alpha_ow.*E_ow.*(T_chg-T_ro))./(1-nu_ow)).*(t_ins./(t_ins+t_out));
% Circumferential stress inner wall [(1)B4,p121]  ()
stressTiw = ((alpha_iw.*E_iw.*(T_chg-T_ro))./(1-nu_iw)).*(t_out./(t_ins+t_out));
% Shear Stress max [(1)B4,p121] ()
stressTs = ( (alpha_avg.*E_avg.*(T_chg-T_ro))./(5.*w_rib.*(1-nu_avg)) ) .* ((t_out.*t_ins)./(t_ins+t_out).^2) .* (w_ch+w_rib);
% Compression Stress max [(1)B4,p121]  ()
stressTc = ( (alpha_avg.*E_avg.*(T_chg-T_ro)) ./ ((w_rib.*r_ow./(w_ch+w_rib)).*(1-nu_avg)) ) .* ((t_out.*t_ins)./(t_ins+t_out));

% PRESSURE STRESS [Pa]
% Maximum bending stress  [(1)A3,p111] ()
stressPb = (deltaP./2).*(w_ch./t_ins).^2;
% Shear stress, ribs (ignore for von mises)
stressPs = (deltaP./2).*(w_ch./t_ins);
% circumferential stress of the engine wall [(1)A4,p112]
stressP_hoop = deltaP.*r_bar./(t_ins+t_out);
% Axial 
stressPa = P_chamber.*A_cc./A_wall;

% {PUT IF STATEMENT HERE}
% if
% If Chamber pressure greater than coolant preasure: Compressive loading, ribs  (ignore for von mises)
% stressPc = (deltaP/w).*(h+w);
% (Expecting this) If Chamber pressure less than coolant preasure: Tensile stress, ribs (ignore for von mises)
stressPt = (deltaP./w_rib).*w_ch;
% end


% Von Mises maximum distortion energy theory
% failure by yielding occurs when, at any point in the body, the distortion energy per unit volume 
% in a state of combined stress becomes equal to that associated with yielding
% Failure occurs when the von Mises stress exceeds the material's 0.2% yield strength.
% https://engrapps.com/mechanical-systems-and-materials/mechanical-components/pressure-vessels/burst-collapse-analysis.php
stressTotalAxial = stressPa;
stressTotalRadial = stressTc + stressPb;
stressTotalHoop = stressTow + stressTiw + stressP_hoop;
stressTotalTow = stressTs+stressPs;

vonMises = sqrt( (1/2).*((stressTotalAxial-stressTotalRadial).^2 + (stressTotalRadial-stressTotalHoop).^2 + (stressTotalHoop-stressTotalAxial).^2 )+3.*stressTotalTow.^2);

end

