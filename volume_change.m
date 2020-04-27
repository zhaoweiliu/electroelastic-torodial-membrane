function [v_ratio] = volume_change(y)
% Compute the enclosed volume change of the toroidal membrane for a given y.
global gamma t_span
area = deformed_area(y);
half_circle = 1/2 * pi * gamma^2;
v_ratio = (area-half_circle)/half_circle;

function [area_deformed] = deformed_area(y)
area_deformed = 0.;
 for j = 1 : length(t_span)-1
     area_deformed = area_deformed + (y(j,3) + y(j+1,3)) * abs(y(j+1,1) - y(j,1))/2;
 end
end

end
