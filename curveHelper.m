function curveHelper(D, t_ins, h_ch, t_out, pos, dtheta, psi)
    inner_wall_hot = [(D/2)' pos' zeros(length(pos),1)];
    writematrix(inner_wall_hot);
    inner_wall_coolant = [(D/2+(t_ins).*cos(psi))' (pos+(t_ins).*sin(-psi))' zeros(length(pos),1)];
    writematrix(inner_wall_coolant);
    outer_wall_coolant = [(D/2+(t_ins+h_ch).*cos(psi))' (pos+(t_ins+h_ch).*sin(-psi))' zeros(length(pos),1)];
    writematrix(outer_wall_coolant);
    outer_wall_ambient = [(D/2+(t_ins+h_ch+t_out).*cos(psi))' (pos+(t_ins+h_ch+t_out).*sin(-psi))' zeros(length(pos),1)];
    writematrix(outer_wall_ambient)

    theta = zeros(1,length(pos));
    theta_0 = -dtheta(1);
    for i=1:length(pos)
        theta(i) = theta_0 + dtheta(i);
        theta_0 = theta(i);
    end
    x = (D/2+(t_ins).*cos(psi)).*cos(theta);
    z = (D/2+(t_ins).*cos(psi)).*sin(theta);
    rib = [x' (pos+(t_ins).*sin(-psi))' z'];
    writematrix(rib)
end