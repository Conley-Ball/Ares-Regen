function [T_c,P_cool,q] = heatTransfer2D(Pc,M,P,A,D,D_h,l_ch,w_ch,num,step,pos,gamma,MW_g,mu_g,Cp_g,C_star,T_aw,k_w,T_inlet,ratio,mdot_f)
    T_c = zeros(1,length(pos));
    P_c = zeros(1,length(pos));
    q = zeros(1,length(pos));

    T_c(1) = T_inlet;
    P_c(1) = Pc;
    
    for i=1:length(pos)
        %% Coolant Properties
        [rho_c,mu_c,k_c,Cp_c,st_c,T_sat_c,v_c,h_lg,rho_c_liquid,rho_c_vapor] = coolant(T_c(i),P_c(i),ratio,mdot_f,l_ch(i),w_ch(i),num);

        %% Coolant Side Heat Transfer

        Re_c = D_h(i)*v_c*rho_c/mu_c;
        Pr_c = mu_c*Cp_c/k_c;
        Nu_c = 0.023*Re_c^0.8*Pr_c^0.4;

        h_c = Nu_c*k_c/D_h(i);



         b(j) = mu_c(i,j) * hf_c(i,j) * (9.81 * (rho_c_lq(i,j) - rho_c_v(i,j)) / st_c(i,j))^0.5 * (cp_c(i,j) / (csf * hf_c(i,j) * Pr_c(i,j)^n_exp))^3;
            
        
        %% Gas Side Heat Transfer
        
        % syms x;
        % a1 = 2*x+3;
        % fun = matlabFunction(a1)
        % ans = fzero(fun,T_c(i))




    end

end