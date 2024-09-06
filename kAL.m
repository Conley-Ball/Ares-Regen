function [k] = kAL(T)
% Material Properties for Additively Manufactured Steel174ph
% Source: Ansys Material Database
%{ There are memory and performance benefits to using "griddedInterpolant" objects over the "interp" functions. 
% griddedInterpolant offers substantial performance improvements for repeated queries of the interpolant object, 
% whereas the interp functions perform a new calculation each time they are called. Also, griddedInterpolant 
% stores the sample points in a memory-efficient format (as a compact grid) and is multithreaded to take 
% advantage of multicore computer processors.
%}

%===PROPERTIES===

% Thermal Conductivity [W/mK]
k_table = [324	179.7;
348	176.2;
375	175.2;
398	175.2;
422	172.7;
449	172.1;
474	171.6;
499	169.7];

% ===Interpolation===
% Thermal Conductivity [W/mK]
k_interp = griddedInterpolant(k_table(:,1), k_table(:,2));
k = k_interp(T);
end
