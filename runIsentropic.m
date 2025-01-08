function [Parameter,Geometry,Gas] = runIsentropic(Parameter,Geometry,Gas)
% RUNISENTROPIC Calculate gas properties from 1-D isentropic flow relations.
% 
%   [Parameter,Geometry,Gas] = runIsentropic(Parameter,Geometry,Gas)
%   calculated flow properties on a set of input parameters in struct
%   Parameters, geometric features in struct Geometry, and gas properties
%   in struct Gas.

% Calculates Mach number using 1-D relations from area ratios using
% flowisentropic().
% Calculation is split into two for loops for supersonic and subsonic
% values.
for lcv=1:Geometry.stations(2)
    Gas.Mach(lcv) = flowisentropic(Gas.specific_heat_ratio(lcv),Geometry.area(lcv)/Geometry.throat_area,'sup');
end
for lcv=Geometry.stations(2)+1:Geometry.stations(4)
    Gas.Mach(lcv) = flowisentropic(Gas.specific_heat_ratio(lcv),Geometry.area(lcv)/Geometry.throat_area,'sub');
end
% Static pressure, total temperature, and static temperature are calculated
% with 1-D isentropic flow relations.
Gas.static_pressure = Parameter.chamber_pressure * (1+(Gas.specific_heat_ratio-1)/2.*Gas.Mach.^2).^(-Gas.specific_heat_ratio/(Gas.specific_heat_ratio-1)); % Pa
Gas.total_temperature = Gas.throat_temperature*Parameter.C_star_efficiency^2*(1+(Gas.specific_heat_ratio-1)/2);
end