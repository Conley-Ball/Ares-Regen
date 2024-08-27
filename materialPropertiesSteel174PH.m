function [E,nu,alpha,k,yield] = materialPropertiesSteel174PH(T)
% Material Properties for Additively Manufactured Steel174ph
% Source: Ansys Material Database
%{ There are memory and performance benefits to using "griddedInterpolant" objects over the "interp" functions. 
% griddedInterpolant offers substantial performance improvements for repeated queries of the interpolant object, 
% whereas the interp functions perform a new calculation each time they are called. Also, griddedInterpolant 
% stores the sample points in a memory-efficient format (as a compact grid) and is multithreaded to take 
% advantage of multicore computer processors.
%}
    

    T = eval(T);
    %===PROPERTIES===
    
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
    
    % Yield Strength [Pa]
    Yield_table = [298	251000000;
    373	232000000;
    423	221000000;
    473	197000000;
    523	148000000;
    1000 138050000];
    
    % ===Interpolation===
    
    % Elastic modulous [Pa]
    E_interp = griddedInterpolant(E_table(:,1), E_table(:,2));
    E = E_interp(T);
    
    % Poison's Ratio
    nu_interp = griddedInterpolant(nu_table(:,1), nu_table(:,2));
    nu = nu_interp(T);
    
    % Thermal Expansion [1/K]
    alpha_interp = griddedInterpolant(alpha_table(:,1), alpha_table(:,2));
    alpha = alpha_interp(T);
    
    % Thermal Conductivity [W/mK]
    k_interp = griddedInterpolant(k_table(:,1), k_table(:,2));
    k = k_interp(T);
    
    % Yield Strength [Pa]
    Yield_interp = griddedInterpolant(Yield_table(:,1), Yield_table(:,2));
    yield = Yield_interp(T);

end