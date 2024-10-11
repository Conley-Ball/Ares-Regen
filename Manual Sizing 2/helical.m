function [x_coords,y_coords,phi] = helical(D,pos,pitch,t_ins,h_ch)
    x_coords = zeros(1,length(pos));
    y_coords = zeros(1,length(pos));
    dr_dz = zeros(1,length(pos));
    theta = zeros(1,length(pos));
    theta_0 = zeros(1,length(pos)); %array intialization
    pos = [pos pos(end)+pos(end)-pos(end-1)]; %adding interpolated value for derivative
    r = D./2+t_ins;
    r = [r r(end)]; %adding interpolated value for derivative
    theta_0(1) = 0;
    alpha = pitch*pi/180; %convert to rad
    for i = 1:length(pos)-1
        dr_dz(i) = (r(i+1)-r(i))/(pos(i+1)-pos(i)); %manually calculate differential
    end
    for i = 1:length(pos)-1
        phi(i) = atan( ((((pos(i+1)-pos(i)) *tan(alpha) )^2 + h_ch(i)^2) / (pos(i+1)-pos(i))^2)^0.5);
        
        theta(i) = theta_0(i) + tan(alpha)/r(i)*sqrt(1+(dr_dz(i))^2)*(pos(i+1)-pos(i)); %angle
        x_coords(i) = r(i)*cos(theta(i));
        y_coords(i) = r(i)*sin(theta(i));
        theta_0(i+1) = atan(y_coords(i)/x_coords(i)); %end point becomes starting point for next point
        if x_coords(i) < 0 %applying arctan to entire domain
            if theta_0(i+1) > -pi
                theta_0(i+1) = theta_0(i+1) - pi;
            end
        end
    end
end

