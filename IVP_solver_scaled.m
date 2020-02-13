function [drhoatpi, etaatpi, tout, y] = IVP_solver_scaled(x)
global epsilon alpha gamma t_span scale_coeff rho0
P = x(1);
y0 = [rho0/scale_coeff,0,0,x(2)];
options = odeset('stats','off','Mass',@mass_matrix,'RelTol',1e-6,'AbsTol',1e-9);
[tout, y] = ode45(@ode_fun,t_span,y0,options);
n = size(y,1);
m = size(tout,1);
if n~= m
    disp('Check. lengths of t and y are not same here!');
end
y_s2 = scale_coeff*y;
drhoatpi = y_s2(n,2);
etaatpi = y_s2(n,3);
% Define mass matrix:

    function Mass = mass_matrix(t,y)
        y_s = scale_coeff*y;
        R = 1 + gamma*cos(t);
        lambda1 =  1/gamma * sqrt(y_s(2)^2 + y_s(4)^2);
        lambda2 =  y_s(1)/R;
        
        A22 = R* (2/(gamma^2) * (1+alpha * lambda2^2)*(1-1/(lambda1^4*lambda2^2))...
            + (8 * y_s(2)^2)/(gamma^4 * lambda1^6 * lambda2^2) ...
            * (1+alpha * lambda2^2) - epsilon * lambda2^2 /(2*gamma^2));
        %
        A24 = R*(8 * y_s(2)*y_s(4))/(gamma^4 * lambda1^6 * lambda2^2) * (1+alpha * lambda2^2);
        
        A42 = A24;
        %
        A44 = R* (2/(gamma^2) * (1+alpha * lambda2^2)*(1-1/(lambda1^4*lambda2^2))...
            + (8 * y_s(4)^2)/(gamma^4 * lambda1^6 * lambda2^2) ...
            * (1+alpha * lambda2^2) - epsilon * lambda2^2 /(2*gamma^2));
        
        Mass = [scale_coeff, 0, 0, 0;
            0, A22*scale_coeff, 0, A24*scale_coeff;
            0, 0, scale_coeff, 0;
            0, A42*scale_coeff, 0, A44*scale_coeff];
    end

% Define ODE:
    function dydt = ode_fun(t,y)
        y_s = scale_coeff*y;
        R = 1 + gamma*cos(t);
        lambda1 =  1/gamma * sqrt(y_s(2)^2 + y_s(4)^2);
        lambda2 =  y_s(1)/R;
        dlambda2dt = y_s(2)/R + (y_s(1)*gamma*sin(t))/(R^2);
        E2 = -(-gamma*sin(t) * ( (2*y_s(2))/(gamma^2) *...
            (1+alpha*lambda2^2) *(1-1/(lambda1^4*lambda2^2)) - ...
            (epsilon * y_s(2) * lambda2^2)/(2*gamma^2)) +...
            R*((2*y_s(2)/(gamma^2)) * (2 * alpha *lambda2 * dlambda2dt) *...
            (1 - 1/(lambda1^4 * lambda2^2)) + ...
            (4*y_s(2) * dlambda2dt)/(gamma^2 * lambda1^4 * lambda2^3)*(1 + alpha*lambda2^2) -...
            (2*epsilon*y_s(2))/(gamma^2) * lambda2 *dlambda2dt) - ...
            2*lambda2*(1 + alpha*lambda1^2)*(1-1/(lambda1^2*lambda2^4)) +...
            (epsilon * lambda2 * lambda1^2)/2 + P*y_s(1)*y_s(4)/gamma);
        %
        E4 = -(-gamma*sin(t) * ((2*y_s(4))/(gamma^2) *...
            (1+alpha*lambda2^2) *(1-1/(lambda1^4*lambda2^2)) - ...
            (epsilon * y_s(4) * lambda2^2)/(2*gamma^2)) +...
            R*((2*y_s(4)/(gamma^2)) * ( 2 * alpha *lambda2 * dlambda2dt) *...
            (1 - 1/(lambda1^4 * lambda2^2)) + ...
            (4*y_s(4) * dlambda2dt)/(gamma^2 * lambda1^4 * lambda2^3)*(1 + alpha*lambda2^2) -...
            (2*epsilon*y_s(4))/(gamma^2) * lambda2 *dlambda2dt) - ...
            (P*y_s(1)*y_s(2) / gamma));
        
        dydt = [y_s(2); E2; y_s(4); E4];
    end

end