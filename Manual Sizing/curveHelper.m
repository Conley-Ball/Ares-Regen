function curveHelper(D, t_ins, h_ch, t_out, pos)
    inner_wall_hot = cat(2,permute(D/2,[2 1]),permute(pos,[2 1]),zeros(numel(pos),1));
    writematrix(inner_wall_hot);
    inner_wall_coolant = cat(2,permute(D/2+t_ins,[2 1]),permute(pos,[2 1]),zeros(numel(pos),1));
    writematrix(inner_wall_coolant);
    outer_wall_coolant = cat(2,permute(D/2+t_ins+h_ch,[2 1]),permute(pos,[2 1]),zeros(numel(pos),1));
    writematrix(outer_wall_coolant);
    outer_wall_ambient = cat(2,permute(D/2+t_ins+h_ch+t_out,[2 1]),permute(pos,[2 1]),zeros(numel(pos),1));
    writematrix(outer_wall_ambient)
end