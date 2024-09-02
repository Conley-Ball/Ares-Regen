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
k_table = [295	110;
323	111;
373	112;
423	114;
473	113;
523	109;
573	116;
673	116;
723	115;
773	109;
803	109;
843	101;
893	54;
913	48;
973	51];

% ===Interpolation===
% Thermal Conductivity [W/mK]
k_interp = griddedInterpolant(k_table(:,1), k_table(:,2));
k = k_interp(T);
end
