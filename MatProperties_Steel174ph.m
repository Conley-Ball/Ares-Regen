function [row,E,E_iw,E_ow,E_avg,nu,nu_iw,nu_ow,nu_avg,alpha,alpha_iw,alpha_ow,alpha_avg,k,Cp,Yield] = MatProperties_Steel174ph(T_iw,T_ow,T)
% Material Properties for Additively Manufactured Steel174ph
% Source: Ansys Material Database
%{ There are memory and performance benefits to using "griddedInterpolant" objects over the "interp" functions. 
% griddedInterpolant offers substantial performance improvements for repeated queries of the interpolant object, 
% whereas the interp functions perform a new calculation each time they are called. Also, griddedInterpolant 
% stores the sample points in a memory-efficient format (as a compact grid) and is multithreaded to take 
% advantage of multicore computer processors.
%}

%===PROPERTIES===
% Density [kg/m^3] (a
row_table = [295	7790];

% Elastic modulous [Pa]
E_table = [298	204000000000;
494	    195000000000;
580	    187000000000;
650	    182000000000;
728	    176000000000;
798	    168000000000;
885	    153000000000;
957	    142000000000;
1067	129000000000;
1162	117000000000];

% Poison's Ratio
nu_table = [298	0.291;
494	    0.295;
580	    0.296;
650	    0.305;
728	    0.316;
798	    0.309;
885	    0.322;
957	    0.332;
1067	0.348;
1162	0.361];

% Thermal Expansion [1/K]
alpha_table = [473	1.1E-05;
573	1.14E-05;
673	1.18E-05;
773	1.2E-05];

% Thermal Conductivity [W/mK]
k_table = [294	15.2;
461	18;
627	19.9;
794	20.6;
961	28.1];

% Specific heat [J/kg]
Cp_table = [350	475;
375	488.4;
400	497.3;
425	505.4;
450	514.7;
460	522;
475	522.8;
500	533.5;
525	547.7;
550	554.8;
600	563.6;
625	572.6;
650	580.5;
675	588.5;
700	597.8;
725	608;
750	634.3;
775	659.4;
795	681.2;
800	687.5;
825	724.6;
850	740.9;
875	757.6;
900	767.2;
925	778.9];

% Yield Strength [Pa]
Yield_table = [298	251000000;
373	232000000;
423	221000000;
473	197000000;
523	148000000];



% ===Interpolation===
% Density [kg/m^3]
row = row_table(1,2);

% Elastic modulous [Pa]
E_interp = griddedInterpolant(E_table(:,1), E_table(:,2));
E = E_interp(T);
E_iw = E_interp(T_iw);
E_ow = E_interp(T_ow);
E_avg = (E_iw+E_ow)/2;

% Poison's Ratio
nu_interp = griddedInterpolant(nu_table(:,1), nu_table(:,2));
nu = nu_interp(T);
nu_iw = nu_interp(T_iw);
nu_ow = nu_interp(T_ow);
nu_avg = (nu_iw+nu_ow)/2;

% Thermal Expansion [1/K]
alpha_interp = griddedInterpolant(alpha_table(:,1), alpha_table(:,2));
alpha = alpha_interp(T);
alpha_iw = alpha_interp(T_iw);
alpha_ow = alpha_interp(T_ow);
alpha_avg = (alpha_iw+alpha_ow)/2;

% Thermal Conductivity [W/mK]
k_interp = griddedInterpolant(k_table(:,1), k_table(:,2));
k = k_interp(T);

% Specific heat [J/kg]
Cp_interp = griddedInterpolant(Cp_table(:,1), Cp_table(:,2));
Cp = Cp_interp(T);

% Yield Strength [Pa]
Yield_interp = griddedInterpolant(Yield_table(:,1), Yield_table(:,2));
Yield = Yield_interp(T);
end
