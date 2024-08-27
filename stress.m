function [stressTotal,stressTow, stressTiw, stressTs, stressTc,stressPb, stressPs, stressPt, stressP_hoopow,stressP_hoopiw, stressPa,vonMises] = stress(T_iw,T_ow,E_iw,E_ow,alpha_iw,alpha_ow,nu_iw,nu_ow,t_iw,t_ow,h,w,d,r_ow,r_iw,deltaP,P_chamber,A_wall,A_cc)
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


% THERMAL STRESS [Pa]
% Circumferential stress outer wall  [(1)B3,p120] ()
stressTow = ((alpha_ow.*E_ow.*(T_iw-T_ow))./(1-nu_ow)).*(t_iw./(t_iw+t_ow));
% Circumferential stress inner wall [(1)B4,p121]  ()
stressTiw = ((alpha_iw.*E_iw.*(T_iw-T_ow))./(1-nu_iw)).*(t_ow/(t_iw+t_ow));
% Shear Stress max [(1)B4,p121] ()
stressTs = ( (alpha_avg.*E_avg.*(T_iw-T_ow))./(5.*w.*(1-nu_avg)) ) .* ((t_ow.*t_iw)/(t_iw+t_ow).^2) .* (h+w);
% Compression Stress max [(1)B4,p121]  ()
stressTc = ( (alpha_avg.*E_avg.*(T_iw-T_ow)) ./ ((w.*r_ow/(h+w)).*(1-nu_avg)) ) .* ((t_ow.*t_iw)/(t_iw+t_ow));

% PRESSURE STRESS: INNER WALL [Pa]
% Maximum bending stress  [(1)A3,p111] ()
stressPb = (deltaP/2).*(h/t_iw)^2;

    % If Chamber pressure greater than coolant preasure: Compressive loading, ribs  (ignore for von mises)
    % stressPc = (deltaP/w).*(h+w);
    % (Expecting this) If Chamber pressure less than coolant preasure: Tensile stress, ribs (ignore for von mises)
    stressPt = (deltaP/w).*h;

% Shear stress, ribs (ignore for von mises)
stressPs = (deltaP/2).*(h/t_iw);
% circumferential stress in the inner and outer walls ()
stressP_hoopow = deltaP.*r_ow./(t_iw+t_ow);
stressP_hoopiw = deltaP.*r_iw./(t_iw+t_ow);
% Axial 
stressPa = P_chamber.*A_cc./A_wall;


% TOTAL STRESS
% these stresses can be simply added together to determine the total stress on an object because of the superposition principle.
stressTotal = stressTow + stressTiw + stressTs + stressTc + stressPb + stressPs + stressP_hoopow + stressP_hoopiw + stressPa;

% Von Mises maximum distortion energy theory
% failure by yielding occurs when, at any point in the body, the distortion energy per unit volume 
% in a state of combined stress becomes equal to that associated with yielding
% Failure occurs when the von Mises stress exceeds the material's 0.2% yield strength.
% https://engrapps.com/mechanical-systems-and-materials/mechanical-components/pressure-vessels/burst-collapse-analysis.php
stressTotalAxial = stressPa;
stressTotalRadial = stressTc + stressPb;
stressTotalHoop = stressTow + stressTiw + stressP_hoopow + stressP_hoopiw;
stressTotalTow = stressTs;

vonMises = sqrt( (1/2).*((stressTotalAxial-stressTotalRadial).^2 + (stressTotalRadial-stressTotalHoop).^2 + (stressTotalHoop-stressTotalAxial).^2 )+3.*stressTotalTow.^2);

end

