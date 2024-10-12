function curveHelper(D, t_ins, h_ch, t_out, pos, dtheta, psi, w_rib, w_ch, num_ch)
    inner_wall_hot = [(D/2)' pos' zeros(length(pos),1)];
    writematrix(inner_wall_hot);
    inner_wall_coolant = [(D/2+(t_ins).*cos(psi))' (pos+(t_ins).*sin(-psi))' zeros(length(pos),1)];
    writematrix(inner_wall_coolant);
    outer_wall_coolant = [(D/2+(t_ins+h_ch).*cos(psi))' (pos+(t_ins+h_ch).*sin(-psi))' zeros(length(pos),1)];
    writematrix(outer_wall_coolant);
    outer_wall_ambient = [(D/2+(t_ins+h_ch+t_out).*cos(psi))' (pos+(t_ins+h_ch+t_out).*sin(-psi))' zeros(length(pos),1)];
    writematrix(outer_wall_ambient)

    theta = zeros(1,length(pos));
    width_theta = zeros(1,length(pos));
    theta_0 = -dtheta(1);
    for i=1:length(pos)
        theta(i) = theta_0 + dtheta(i);
        theta_0 = theta(i);
        width_theta(i) = pi/num_ch*w_rib(i)/(w_rib(i)+w_ch(i));
    end
    inj_add = 0.0088138;
    helix_1_x_add = (D(end)/2+(t_ins(end)-0.0001).*cos(psi(end))).*cos(theta(end)-width_theta(end));
    helix_1_y_end = (D(end)/2+(t_ins(end)-0.0001).*cos(psi(end))).*sin(theta(end)-width_theta(end));
    helix_2_x_add = (D(end)/2+(t_ins(end)-0.0001).*cos(psi(end))).*cos(theta(end)+width_theta(end));
    helix_2_y_end = (D(end)/2+(t_ins(end)-0.0001).*cos(psi(end))).*sin(theta(end)+width_theta(end));
    helix_3_x_add = (D(end)/2+(t_ins(end)+h_ch(end)+0.0001).*cos(psi(end))).*cos(theta(end)-width_theta(end));
    helix_3_y_end = (D(end)/2+(t_ins(end)+h_ch(end)+0.0001).*cos(psi(end))).*sin(theta(end)-width_theta(end));
    helix_4_x_add = (D(end)/2+(t_ins(end)+h_ch(end)+0.0001).*cos(psi(end))).*cos(theta(end)+width_theta(end));
    helix_4_y_end = (D(end)/2+(t_ins(end)+h_ch(end)+0.0001).*cos(psi(end))).*sin(theta(end)+width_theta(end));
    pos_add_12 = (pos(end)+(t_ins(end)-0.0001).*sin(-psi(end))) + inj_add;
    pos_add_34 = (pos(end)+(t_ins(end)+h_ch(end)+0.0001).*sin(-psi(end))) + inj_add;
    helix_1 = [[(D/2+(t_ins-0.0001).*cos(psi)).*cos(theta-width_theta) helix_1_x_add]' [(pos+(t_ins-0.0001).*sin(-psi)) pos_add_12]' [(D/2+(t_ins-0.0001).*cos(psi)).*sin(theta-width_theta) helix_1_y_end]'];
    helix_2 = [[(D/2+(t_ins-0.0001).*cos(psi)).*cos(theta+width_theta) helix_2_x_add]' [(pos+(t_ins-0.0001).*sin(-psi)) pos_add_12]' [(D/2+(t_ins-0.0001).*cos(psi)).*sin(theta+width_theta) helix_2_y_end]'];
    helix_3 = [[(D/2+(t_ins+h_ch+0.0001).*cos(psi)).*cos(theta-width_theta) helix_3_x_add]' [(pos+(t_ins+h_ch+0.0001).*sin(-psi)) pos_add_34]' [(D/2+(t_ins+h_ch+0.0001).*cos(psi)).*sin(theta-width_theta) helix_3_y_end]'];
    helix_4 = [[(D/2+(t_ins+h_ch+0.0001).*cos(psi)).*cos(theta+width_theta) helix_4_x_add]' [(pos+(t_ins+h_ch+0.0001).*sin(-psi)) pos_add_34]' [(D/2+(t_ins+h_ch+0.0001).*cos(psi)).*sin(theta+width_theta) helix_4_y_end]'];
    writematrix(helix_1)
    writematrix(helix_2)
    writematrix(helix_3)
    writematrix(helix_4)
    % x = (D/2+(t_ins).*cos(psi)).*cos(theta);
    % z = (D/2+(t_ins).*cos(psi)).*sin(theta);
    % rib = [x' (pos+(t_ins).*sin(-psi))' z'];
end