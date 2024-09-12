function [E,nu,alpha,k,Cp,Yield] = materialProperties(T,material)
    if strcmp(material,'aluminum') % AlSi10Mg 
        % Elastic modulous [Pa]
        E_table = [298	76600000000;
        323	76100000000;
        373	74300000000;
        423	72700000000;
        473	70600000000;
        523	68900000000;
        573	67000000000];
        
        % Poison's Ratio
        nu_table = [298	0.33;
        323	0.33;
        373	0.33;
        423	0.33;
        473	0.33;
        523	0.33;
        573	0.33;];
        
        % Thermal Expansion [1/K]
        alpha_table = [373	2.06E-05;
        423	2.36E-05;
        473	2.47E-05;
        523	2.58E-05;
        573	3.04E-05;
        623	3.29E-05;
        673	2.71E-05;
        723	2.44E-05];
        
        % Thermal Conductivity [W/mK]
        k_table = [324	179.7;
        348	176.2;
        375	175.2;
        398	175.2;
        422	172.7;
        449	172.1;
        474	171.6;
        499	169.7];
        
        % Specific heat [J/kg]
        Cp_table = [295	91
        423	1010;
        548	1025;
        698	1136;
        823	1326];
        
        % Yield Strength [Pa]
        Yield_table = [298	251000000;
        373	232000000;
        423	221000000;
        473	197000000;
        523	148000000];
    elseif strcmp(material,'steel') % Steel174ph
        % Elastic modulous [Pa]
        E_table = [298	204000000000;
        494	195000000000;
        580	187000000000;
        650	182000000000;
        728	176000000000;
        798	168000000000;
        885	153000000000;
        957	142000000000;
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
        alpha_table = [373	1.04E-05;
        473	1.10E-05;
        573	1.14E-05;
        673	1.18E-05;
        773	1.20E-05];
        
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
        Yield_table = [473	860000000;
        698	754000000;
        748	700000000;
        803	621000000;
        858	508000000;
        918	371000000];

    elseif strcmp(material,'inconel') % Inconel 718
        % Elastic modulous [Pa]
        E_table = [294	165000000000;
        810	    152000000000;
        1088	110000000000;
        1255	55000000000;
        1366	34000000000];
        
        % Poison's Ratio
        nu_table = [294	0.3;
        810	    0.28;
        1088	0.323;
        1255	0.368;
        1366	0.4];
        
        % Thermal Expansion [1/K]
        alpha_table = [700	1.44E-05;
        811	    1.49E-05;
        922	    1.54E-05;
        1144	1.75E-05;
        1366	1.83E-05;
        2100	1.83E-05];
        
        % Thermal Conductivity [W/mK]
        k_table = [295	11.9;
        506	    13.7;
        721	    16.9;
        930	    21.7;
        1139	25.6;
        1352	22.9;
        1562	19.1;
        1773	17.7];
        
        % Specific heat [J/kg]
        Cp_table = [293	421;
        373	    442;
        473	    453;
        573	    472;
        673	    481;
        773	    502;
        873	    527;
        973	    562;
        1073	606;
        1123	628;
        1173	636;
        1273	647;
        1373	651;
        1773	652];
        
        % Yield Strength [Pa]
        Yield_table = [294	648000000;
        811	    558000000;
        1089	338000000;
        1255	90000000;
        1366	34000000];
    else
        error('Invalid material')
    end
    
    % Elastic modulous [Pa]
    E_interp = griddedInterpolant(E_table(:,1), E_table(:,2));
    E = E_interp(T);
    % E = interp1(E_table(:,1), E_table(:,2),T,'pchip');
    
    % Poison's Ratio
    nu_interp = griddedInterpolant(nu_table(:,1), nu_table(:,2));
    nu = nu_interp(T);
    % nu = interp1(nu_table(:,1), nu_table(:,2),T,'pchip');
    
    % Thermal Expansion [1/K]
    alpha_interp = griddedInterpolant(alpha_table(:,1), alpha_table(:,2));
    alpha = alpha_interp(T);
    % alpha = interp1(alpha_table(:,1), alpha_table(:,2),T,'pchip');
    
    % Thermal Conductivity [W/mK]
    k_interp = griddedInterpolant(k_table(:,1), k_table(:,2));
    k = k_interp(T);
    % k = interp1(k_table(:,1), k_table(:,2),T,'pchip');
    
    % Specific heat [J/kg]
    Cp_interp = griddedInterpolant(Cp_table(:,1), Cp_table(:,2));
    Cp = Cp_interp(T);
    % Cp = interp1(Cp_table(:,1), Cp_table(:,2),T,'pchip');
    
    % Yield Strength [Pa]
    Yield_interp = griddedInterpolant(Yield_table(:,1), Yield_table(:,2));
    Yield = Yield_interp(T)*1;
    % Yield = interp1(Yield_table(:,1), Yield_table(:,2),T,'pchip');
end