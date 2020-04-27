function [s11,s22,is_negative] = principle_stresses(y,x,t)
% Compute the principle stresses at all thetas.
global epsilon alpha gamma t_span H R_b;
lambda1 = sqrt(y(:,2).^2 + y(:,4).^2)/gamma;
lambda2 = y(:,1)./(1+gamma*cos(t));
lambda3 = 1./(lambda1.*lambda2);
s11 = - (x(1) * H)/R_b + 2 .* lambda1 .* lambda1 +...
    2*alpha .* lambda1 .* lambda1 .* (lambda2.*lambda2 + 1./(lambda1 .* lambda1 .* lambda2 .* lambda2))...
    -2./(lambda1.*lambda1.*lambda2.*lambda2) - 2*alpha.*(1./(lambda1.*lambda1) + 1./(lambda2.*lambda2))...
    -0.5 * epsilon .* lambda1 .* lambda1 .* lambda2 .* lambda2;
s22 = - (x(1) * H)/R_b + 2 .* lambda2 .* lambda2 +...
    2*alpha .* lambda2 .* lambda2 .* (lambda1.*lambda1 + 1./(lambda1 .* lambda1 .* lambda2 .* lambda2))...
    -2./(lambda1.*lambda1.*lambda2.*lambda2) - 2*alpha.*(1./(lambda1.*lambda1) + 1./(lambda2.*lambda2))...
    -0.5 * epsilon .* lambda1 .* lambda1 .* lambda2 .* lambda2;
is_negative = 0;
v_negative_stress = 0;
p_negative_stress = 0;
for ie = 1:length(s22)
    if s22(ie) < 0
        is_negative = 1;
    end
end
end

