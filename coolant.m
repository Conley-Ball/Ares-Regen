function [rho_c,T_sat_c,mu_c,k_c,Cp_c,v_c] = coolant(T_c,P_c,ratio,mdot_f,w_ch,h_ch,num)   
        rho_water = py.CoolProp.CoolProp.PropsSI('D','T', T_c, 'P', P_c, 'water');
        mu_water = py.CoolProp.CoolProp.PropsSI('V','T', T_c, 'P', P_c, 'water');
        k_water = py.CoolProp.CoolProp.PropsSI('L','T', T_c, 'P', P_c, 'water');
        Cp_water = py.CoolProp.CoolProp.PropsSI('C','T', T_c, 'P', P_c, 'water');
        T_sat_water = py.CoolProp.CoolProp.PropsSI('T', 'P', P_c, 'Q', 0, 'water');

        rho_ethanol = py.CoolProp.CoolProp.PropsSI('D','T', T_c, 'P', P_c, 'ethanol');
        mu_ethanol = py.CoolProp.CoolProp.PropsSI('V','T', T_c, 'P', P_c, 'ethanol');
        k_ethanol = py.CoolProp.CoolProp.PropsSI('L','T', T_c, 'P', P_c, 'ethanol');
        Cp_ethanol = py.CoolProp.CoolProp.PropsSI('C','T', T_c, 'P', P_c, 'ethanol');
        T_sat_ethanol = py.CoolProp.CoolProp.PropsSI('T', 'P', P_c, 'Q', 0, 'ethanol');

        rho_c = ((1-ratio)/rho_water+ratio/rho_ethanol)^-1;
        mu_c = mu_water*(1-ratio)+mu_ethanol*ratio;
        k_c = k_water*(1-ratio)+k_ethanol*ratio;
        Cp_c = Cp_water*(1-ratio)+Cp_ethanol*ratio;
        T_sat_c = T_sat_water*(1-ratio)+T_sat_ethanol*ratio;
        v_c = mdot_f/((w_ch*h_ch)*num*rho_c);


end