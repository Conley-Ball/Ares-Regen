function [Thermal] = Bartz(Parameter,Geometry,Gas,Coolant,Thermal,T_w)
% BARTZ evaluate 2-D conduction and convection using Bartz correlation.
%
%   [Thermal] = Bartz(Parameter,Geometry,Gas,Coolant,Thermal,T_w) finds
%   temperatures and heats of the engine wall by setting the function x to
%   zero.
%
%   See also: RUNHEATTRANSFER

% Index is stored as local variable for convenience.
I = Parameter.index;

%% Coolant

% Coefficient S is found for Chen correlation [3].
S = (1+2.53e-6*(Coolant.Reynolds_number(I))^1.17)^-1;

% If the wall temperature exceeds the saturation temperature, nucleate
% boiling is accounted for
if T_w(3) > Coolant.saturation_temperature(I)

    % Pressure difference is found between saturation pressure at wall
    % temperature and saturation pressure at bulk temperature.
    if T_w(3) < 514
    delta_P_sat = py.CoolProp.CoolProp.PropsSI('P', 'T', real(T_w(3)), 'Q', 1, 'water')*(1-Parameter.ethanol_ratio)+py.CoolProp.CoolProp.PropsSI('P', 'T', real(T_w(3)), 'Q', 1, 'ethanol')*Parameter.ethanol_ratio-Coolant.saturation_pressure(I);
    
    % If the wall temperature is above the critical point, there is no
    % saturation pressure, thus delta P is simply maximized.
    else
    delta_P_sat = py.CoolProp.CoolProp.PropsSI('P', 'T', 514, 'Q', 1, 'water')*(1-Parameter.ethanol_ratio)+py.CoolProp.CoolProp.PropsSI('P', 'T', 514, 'Q', 1, 'ethanol')*Parameter.ethanol_ratio-Coolant.saturation_pressure(I); 
    end

    % Temperature difference between wall and saturation temperature is
    % found.
    delta_T_sat = T_w(3) - Coolant.saturation_temperature(I);

    % Nucleate boiling coefficient is found by multiplying the first term
    % with the requisite delta P and T terms [3].
    Coolant.heat_transfer_coefficient_nucleate_boiling = Coolant.heat_transfer_coefficient_nucleate_boiling*delta_T_sat^0.25*delta_P_sat^0.75;

    % Final coolant convection coefficient is found by combining the single
    % phase and nucleate boiling coefficients [3].
    Thermal.h_c(I) = Coolant.heat_transfer_coefficient_single_phase + S*Coolant.heat_transfer_coefficient_nucleate_boiling*(T_w(3) - Coolant.saturation_temperature(I))/(T_w(3) - Coolant.temperature(I));

% If the wall temperature does not exceed the bulk saturation temperature,
% nucleate boiling does not occur.
else
    Thermal.h_c(I) = Coolant.heat_transfer_coefficient_single_phase;
end

%% Bartz correlation

% Auxiliary variable sigma is found for the inner wall in near the channel
% and rib.
sigma(1) = ( (0.5 * T_w(1)/Gas.total_temperature(I) * (1+(Gas.specific_heat_ratio(I)-1)/2*Gas.Mach(I)^2) + 0.5)^0.68 * (1+(Gas.specific_heat_ratio(I)-1)/2*Gas.Mach(I)^2)^0.12 )^-1;
sigma(2) = ( (0.5 * T_w(2)/Gas.total_temperature(I) * (1+(Gas.specific_heat_ratio(I)-1)/2*Gas.Mach(I)^2) + 0.5)^0.68 * (1+(Gas.specific_heat_ratio(I)-1)/2*Gas.Mach(I)^2)^0.12 )^-1;

% Convective heat transfer coefficient is found for channel and rib
% section.
Thermal.h_g(1,I) = (0.026/Geometry.throat_diameter^0.2 * (Gas.viscosity(I)^0.2*Gas.specific_heat_capacity(I)/Gas.Prandtl_number(I)^0.6) * (Parameter.chamber_pressure/(Parameter.characteristic_velocity*Parameter.C_star_efficiency))^0.8 * (Geometry.throat_diameter/(0.5258*0.5*Geometry.throat_diameter))^0.1 * (Geometry.throat_area/Geometry.area(I))^0.9) * sigma(1);
Thermal.h_g(2,I) = (0.026/Geometry.throat_diameter^0.2 * (Gas.viscosity(I)^0.2*Gas.specific_heat_capacity(I)/Gas.Prandtl_number(I)^0.6) * (Parameter.chamber_pressure/(Parameter.characteristic_velocity*Parameter.C_star_efficiency))^0.8 * (Geometry.throat_diameter/(0.5258*0.5*Geometry.throat_diameter))^0.1 * (Geometry.throat_area/Geometry.area(I))^0.9) * sigma(2);

%% Geometry
%  Geometric dimensions of the wall are stored in local variables in
%  accordance with naming convention in [5].
t = Geometry.inner_wall_thickness(I);
t_r = t/2;
h = Geometry.channel_height(I);
d = Geometry.outer_wall_thickness(I);
d_r = d/2;
w = Geometry.rib_width(I);
b = Geometry.channel_width(I);
l = (b+w)/2;

% Auxiliary variables for fin heat transfer are calculated [5].
m = (2*Thermal.h_c(I)/(Parameter.k_inner_channel_temperature(2)*w))^0.5;
epsilon = (2*Parameter.k_inner_channel_temperature(2)/(Thermal.h_c(I)*w))^0.5;
H = (d_r/Parameter.k_inner_channel_temperature(2) + w*l/(2*d*Parameter.k_inner_channel_temperature(2)) + w*d_r/(b*Parameter.k_inner_channel_temperature(2)) + w/(b*Thermal.h_c(I)))^-1;
beta = (H/(Parameter.k_inner_channel_temperature(2)*m)*sinh(m*h)+cosh(m*h))^-1;
eta = (cosh(m*h)-beta)/sinh(m*h);

% Inner wall temperature corresponding to the location between the coating
% and the metal is found [5].
Thermal.inner_wall_temperature(1,I) = T_w(1) - Parameter.res*Thermal.h_g(1,I)*(Gas.adiabatic_wall_temperature(I)-T_w(1));
Thermal.inner_wall_temperature(2,I) = T_w(2) - Parameter.res*Thermal.h_g(2,I)*(Gas.adiabatic_wall_temperature(I)-T_w(2));

% Temperature in the centre of the inner wall is found [5].
Thermal.inner_temperature(1,I) = Thermal.inner_wall_temperature(1,I) - t_r/Parameter.k_inner_temperature(1)/Parameter.res*(T_w(1)-Thermal.inner_wall_temperature(1,I));
Thermal.inner_temperature(2,I) = Thermal.inner_wall_temperature(2,I) - t_r/Parameter.k_inner_temperature(2)/Parameter.res*(T_w(2)-Thermal.inner_wall_temperature(2,I));

% Temperature in the inner wall adjacenct to the channel and rib is found.
Thermal.inner_channel_temperature(1,I) = (Parameter.k_inner_channel_temperature(1)/(t-t_r)*Thermal.inner_temperature(1,I)+Thermal.h_c(I)*Coolant.temperature(I))/(Thermal.h_c(I)+Parameter.k_inner_channel_temperature(1)/(t-t_r));
Thermal.inner_channel_temperature(2,I) = (Parameter.k_inner_channel_temperature(2)/(t-t_r)*Thermal.inner_temperature(2,I)+eta*epsilon*Thermal.h_c(I)*Coolant.temperature(I))/(eta*epsilon*Thermal.h_c(I)+Parameter.k_inner_channel_temperature(2)/(t-t_r));

% Internal heat fluxes are found and named according to the naming
% convention in [5].
Q_cb = b*Thermal.h_c(I)*(Thermal.inner_channel_temperature(1,I)-Coolant.temperature(I)); % W/m
Q_rb = w*epsilon*eta*Thermal.h_c(I)*(Thermal.inner_channel_temperature(2,I)-Coolant.temperature(I)); % W/m
Q_cr = t*0.5*(Parameter.k_inner_temperature(1)+Parameter.k_inner_temperature(2))/l*(Thermal.inner_temperature(1,I)-Thermal.inner_temperature(2,I)); % W/m
Q_chg = b*Thermal.h_g(1,I)*(Gas.adiabatic_wall_temperature(I)-T_w(1)); % W/m
Q_rhg = w*Thermal.h_g(2,I)*(Gas.adiabatic_wall_temperature(I)-T_w(2)); % W/m

% Total heat flux (per unit area) is found from the heat fluxes.
Thermal.heat_flux(I) = (Q_chg+Q_rhg)/(b+w);

% The variables x will be minimized by fsolve.

% Indices 1 and 2 represent conservation across the two control volumes.
Thermal.x(1) = -Q_chg+2*Q_cr+Q_cb;
Thermal.x(2) = -Q_rhg-2*Q_cr+Q_rb;
% Index 3 represents that the inner channel temperature is equivalent to
% the temeprature used earlier for nucleate boiling calculations.
Thermal.x(3) = Thermal.inner_channel_temperature(1,I) - T_w(3);

end

% References:
% [3] Fang, Yuan, et al. Review of correlations for subcooled flow
% [5] Fagherazzi, Matteo et al. A Simplified Thermal Analysis Model for Regeneratively Cooled Rocket Engine Thrust Chambers