function [rho_c,T_sat_c,P_sat_c,mu_c,k_c,Cp_c,st_c,v_c,D_h,rho_c_liquid,rho_c_vapor,h_c_liquid,h_c_vapor] = coolant(T_c,P_c,ratio,mdot_f,w_ch,h_ch,num,fillet,alpha)   
        rho_water = py.CoolProp.CoolProp.PropsSI('D','T', T_c, 'P', P_c, 'water');
        mu_water = py.CoolProp.CoolProp.PropsSI('V','T', T_c, 'P', P_c, 'water');
        k_water = py.CoolProp.CoolProp.PropsSI('L','T', T_c, 'P', P_c, 'water');
        Cp_water = py.CoolProp.CoolProp.PropsSI('C','T', T_c, 'P', P_c, 'water');
        st_water = py.CoolProp.CoolProp.PropsSI('I','P', P_c,'Q',0, 'water');
        T_sat_water = py.CoolProp.CoolProp.PropsSI('T', 'P', P_c, 'Q', 0, 'water');
        P_sat_water = py.CoolProp.CoolProp.PropsSI('P', 'T', T_c, 'Q', 0, 'water');

        rho_ethanol = py.CoolProp.CoolProp.PropsSI('D','T', T_c, 'P', P_c, 'ethanol');
        mu_ethanol = py.CoolProp.CoolProp.PropsSI('V','T', T_c, 'P', P_c, 'ethanol');
        k_ethanol = py.CoolProp.CoolProp.PropsSI('L','T', T_c, 'P', P_c, 'ethanol');
        Cp_ethanol = py.CoolProp.CoolProp.PropsSI('C','T', T_c, 'P', P_c, 'ethanol');
        st_ethanol = py.CoolProp.CoolProp.PropsSI('I','P', P_c,'Q',0, 'ethanol');
        T_sat_ethanol = py.CoolProp.CoolProp.PropsSI('T', 'P', P_c, 'Q', 0, 'ethanol');
        P_sat_ethanol = py.CoolProp.CoolProp.PropsSI('P', 'T', T_c, 'Q', 0, 'ethanol');

        rho_water_liquid = py.CoolProp.CoolProp.PropsSI('D', 'P', P_c, 'Q', 0, 'water');
        h_water_liquid = py.CoolProp.CoolProp.PropsSI('H', 'P', P_c, 'Q', 0, 'water');
        rho_water_vapor = py.CoolProp.CoolProp.PropsSI('D', 'P', P_c, 'Q', 1, 'water');
        h_water_vapor = py.CoolProp.CoolProp.PropsSI('H', 'P', P_c, 'Q', 1, 'water');
        
        rho_ethanol_liquid = py.CoolProp.CoolProp.PropsSI('D', 'P', P_c, 'Q', 0, 'ethanol');
        h_ethanol_liquid = py.CoolProp.CoolProp.PropsSI('H', 'P', P_c, 'Q', 0, 'ethanol');
        rho_ethanol_vapor = py.CoolProp.CoolProp.PropsSI('D', 'P', P_c, 'Q', 1, 'ethanol');
        h_ethanol_vapor = py.CoolProp.CoolProp.PropsSI('H', 'P', P_c, 'Q', 1, 'ethanol');
        
        rho_c_liquid = ((1-ratio)/rho_water_liquid+ratio/rho_ethanol_liquid)^-1;
        rho_c_vapor = ((1-ratio)/rho_water_vapor+ratio/rho_ethanol_vapor)^-1;
        h_c_liquid = h_water_liquid*(1-ratio)+h_ethanol_liquid*ratio;
        h_c_vapor = h_water_vapor*(1-ratio)+h_ethanol_vapor*ratio;

        rho_c = ((1-ratio)/rho_water+ratio/rho_ethanol)^-1;
        mu_c = mu_water*(1-ratio)+mu_ethanol*ratio;
        k_c = k_water*(1-ratio)+k_ethanol*ratio;
        Cp_c = Cp_water*(1-ratio)+Cp_ethanol*ratio;
        st_c = st_water*(1-ratio)+st_ethanol*ratio;
        T_sat_c = T_sat_water*(1-ratio)+T_sat_ethanol*ratio;
        P_sat_c = P_sat_water*(1-ratio)+P_sat_ethanol*ratio;

        A_c = w_ch*h_ch - (2*fillet)^2 + pi*fillet^2;
        Per = 2*(w_ch-2*fillet)+2*(h_ch-2*fillet)+2*fillet*pi;
        D_h = 4*A_c/Per;
        v_c = mdot_f/(A_c*num*rho_c*cos(alpha));
end