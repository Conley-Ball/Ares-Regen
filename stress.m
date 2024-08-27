function [stress_total,stress_T_ow, stress_T_iw, stress_T_s, stress_T_c,stress_P_b, stress_P_c, stress_P_s, stress_P_t, stress_P_hoop, stress_P_a,yield_out] = stress(T_ci,T_co,P,P_c,A,D,b,w,t,h,d,num_ch)


    t = t*ones(1,length(A));
    d = d*ones(1,length(A));
    w = w*ones(1,length(A));
    h = h*ones(1,length(A));
    for i=1:length(A)
        % MATREIAL PROPERTIES
        [E_iw,nu_iw,alpha_iw,k_iw,yield_iw] = materialPropertiesSteel174PH(T_ci(i));
        [E_ow,nu_ow,alpha_ow,k_ow,yield_ow] = materialPropertiesSteel174PH(T_co(i));
        E = 0.5*(E_iw+E_ow);
        nu = 0.5*(nu_iw+nu_ow);
        alpha = 0.5*(alpha_iw+alpha_ow);
        k = 0.5*(k_iw+k_ow);
        yield = 0.5*(yield_iw+yield_ow);
        
        % THERMAL STRESS [Pa]
        % Circumferential stress outer wall  [(1)B3,p120] 
        stress_T_ow(i) = ((alpha_ow*E_ow*(T_ci(i)-T_co(i)))/(1-nu_ow))*(t(i)/(t(i)+d(i)));
        % Circumferential stress inner wall [(1)B4,p121] 
        stress_T_iw(i) = ((alpha_iw*E_iw*(T_ci(i)-T_co(i)))/(1-nu_iw))*(d(i)/(t(i)+d(i)));
        % Shear Stress max [(1)B4,p121]
        stress_T_s(i) = ( (alpha*E*(T_ci(i)-T_co(i)))/(5*w(i)*(1-nu)) ) * ((t(i)*d(i))/(t(i)+d(i))^2) * (b(i)+w(i));
        % Compression Stress max [(1)B4,p121]
        stress_T_c(i) = ( (alpha*E*(T_ci(i)-T_co(i))) / ((w(i)*(D(i)/2+t(i)+h(i)+d(i))/(b(i)+w(i)))*(1-nu))) * ((d(i)*t(i))/(t(i)+d(i)));
        
        % PRESSURE STRESS: INNER WALL [Pa]
        % Maximum bending stress  [(1)A3,p111]
        stress_P_b(i) = ((P_c(i)-P(i))/2)*(b(i)/t(i))^2;
        % Compressive loading, ribs
        stress_P_c(i) = ((P_c(i)-P(i))/w(i))*(b(i)+w(i));
        % Tensile stress, ribs
        stress_P_t(i) = ((P_c(i)-P(i))/w(i))*b(i);
        % Shear stress, ribs
        stress_P_s(i) = ((P_c(i)-P(i))/2)*(b(i)/t(i));
        % circumferential stress in the inner and outer walls
        stress_P_hoop(i) = (P_c(i)-P(i))*(D(i)/2+t(i)+h(i)+d(i))/(t(i)+d(i));
        % Hoop 
        stress_P_hoop(i) = (P_c(i)+P(i))*(D(i)/2)/(t(i)+d(i));
        % Axial 
        stress_P_a(i) = (P_c(i)-P(i))*A(i)/(pi*(D(i)/2+t(i)+h(i)+d(i))^2-A(i)-(num_ch*b(i)*h(i)));
        % if stress_P_a(i) < 0
        %     stress_P_a(i) = 0;
        % end
        yield_out(i) = yield;
    
    end
    
    
    % TOTAL STRESS
    % these stresses can be simply added together to determine the total stress on an object because of the superposition principle.
    stress_total = stress_T_ow + stress_T_iw + stress_T_s + stress_T_c + stress_P_b + stress_P_c + stress_P_s + stress_P_t + stress_P_hoop + stress_P_a;
    
    % Von Mises maximum distortion energy theory
    % failure by yielding occurs when, at any point in the body, the distortion energy per unit volume 
    % in a state of combined stress becomes equal to that associated with yielding

end

