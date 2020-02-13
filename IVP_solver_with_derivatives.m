function [drhoatpi, etaatpi, tout, y, dy] = IVP_solver_with_derivatives(x)
global epsilon alpha gamma t_span scale_coeff rho0
P = x(1);
y0 = [rho0;0;0;x(2)*scale_coeff];
options = odeset('stats','off','Mass',@mass_matrix,'RelTol',1e-4,'AbsTol',1e-6);
[tout, y] = ode45(@ode_fun,t_span,y0,options);
n = size(y,1);
m = size(tout,1);
if n~= m
    disp('Check. lengths of t and y are not same here!');
end
drhoatpi = y(n,2);
etaatpi = y(n,3);
dy = zeros(length(t_span), 4);
for i = 1:length(t_span)
    mass_sol = mass_matrix(t_span(i),y(i,:));
    y_sol = ode_fun(t_span(i),y(i,:));
    dy(i,:) = (mass_sol \ y_sol)';
end

% Define mass matrix
    function Mass = mass_matrix(t,y)
        
        R = 1 + gamma*cos(t);      
        lambda1 = 1/gamma * sqrt(y(2)^2 + y(4)^2);
        lambda2 = y(1)/R;  
        
        A22 = R* (2/(gamma^2) * (1+alpha * lambda2^2)*(1-1/(lambda1^4*lambda2^2))...
            + (8 * y(2)^2)/(gamma^4 * lambda1^6 * lambda2^2) ...
            * (1+alpha * lambda2^2) - epsilon * lambda2^2 /(2*gamma^2));         
        A24 = R*(8 * y(2)*y(4))/(gamma^4 * lambda1^6 * lambda2^2) * (1+alpha * lambda2^2);
        
        A42 = A24;         
        A44 = R* (2/(gamma^2) * (1+alpha * lambda2^2)*(1-1/(lambda1^4*lambda2^2))...
            + (8 * y(4)^2)/(gamma^4 * lambda1^6 * lambda2^2) ...
            * (1+alpha * lambda2^2) - epsilon * lambda2^2 /(2*gamma^2));
        
        Mass = [1, 0, 0, 0;
                0, A22, 0, A24;
                0, 0, 1, 0;
                0, A42, 0, A44];
    end

% Define ODE:
    function dydt = ode_fun(t,y)
        R = 1 + gamma*cos(t);
        lambda1 = 1/gamma * sqrt(y(2)^2 + y(4)^2);
        lambda2 = y(1)/R;
        dlambda2dt = y(2)/R + (y(1)*gamma*sin(t))/(R^2);
                 
        E2 = -( -gamma*sin(t) * ( (2*y(2))/(gamma^2) *...
            (1+alpha*lambda2^2) *(1-1/(lambda1^4*lambda2^2)) - ...
            (epsilon * y(2) * lambda2^2)/(2*gamma^2)) +...
            R*((2*y(2)/(gamma^2)) * (2 * alpha *lambda2 * dlambda2dt) *...
            (1 - 1/(lambda1^4 * lambda2^2)) + ...
            (4*y(2) * dlambda2dt)/(gamma^2 * lambda1^4 * lambda2^3)*(1 + alpha*lambda2^2) -...
            (2*epsilon*y(2))/(gamma^2) * lambda2 *dlambda2dt) - ...
            2*lambda2*(1 + alpha*lambda1^2)*(1-1/(lambda1^2*lambda2^4)) +...
            (epsilon * lambda2 * lambda1^2)/2 + P*y(1)*y(4)/gamma);
                
        E4 = -(-gamma*sin(t) * ((2*y(4))/(gamma^2) *...
            (1+alpha*lambda2^2) *(1-1/(lambda1^4*lambda2^2)) - ...
            (epsilon * y(4) * lambda2^2)/(2*gamma^2)) +...
            R*((2*y(4)/(gamma^2)) * ( 2 * alpha *lambda2 * dlambda2dt) *...
            (1 - 1/(lambda1^4 * lambda2^2)) + ...
            (4*y(4) * dlambda2dt)/(gamma^2 * lambda1^4 * lambda2^3)*(1 + alpha*lambda2^2) -...
            (2*epsilon*y(4))/(gamma^2) * lambda2 *dlambda2dt) - ...
            (P*y(1)*y(2) / gamma));
        
        dydt = [y(2); E2; y(4); E4];
    end
end