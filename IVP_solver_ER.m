function [drhoatpi, etaatpi, tout, y] = IVP_solver_ER(x)
global epsilon alpha gamma t_span rho0 H R_b
P = x(1);
y0 = [rho0,0,0,x(2)];
options = odeset('stats','off','Mass',@mass_matrix,'RelTol',1e-3,'AbsTol',1e-6);
[tout, y] = ode45(@ode_fun,t_span,y0,options);
n = size(y,1);
m = size(tout,1);
if n~= m
    disp('Check. lengths of t and y are not same here!');
end
drhoatpi = y(n,2);
etaatpi = y(n,3);

    function Mass = mass_matrix(t,y)
        lambda1 = 1/gamma * sqrt(y(2)^2 + y(4)^2);
        R = 1 + gamma*cos(t);
        lambda2 = y(1)/R;
        s2 = - (P * H)/R_b + 2 * lambda2 * lambda2 +...
            2*alpha * lambda2 * lambda2 * (lambda1*lambda1 + 1/(lambda1 * lambda1 * lambda2 * lambda2))...
            -2/(lambda1*lambda1*lambda2*lambda2) - 2*alpha*(1/(lambda1*lambda1) + 1/(lambda2*lambda2))...
            -0.5 * epsilon * lambda1 * lambda1 * lambda2 * lambda2;
        if s2 > 0
            Mass = mass_matrix1(t,y);
        else
            Mass = mass_matrix2(t,y);
            
        end
    end

    function dydt = ode_fun(t,y)
        lambda1 = 1/gamma * sqrt(y(2)^2 + y(4)^2);
        R = 1 + gamma*cos(t);
        lambda2 = y(1)/R;
        s2 = - (P * H)/R_b + 2 * lambda2 * lambda2 +...
            2*alpha * lambda2 * lambda2 * (lambda1*lambda1 + 1/(lambda1 * lambda1 * lambda2 * lambda2))...
            -2/(lambda1*lambda1*lambda2*lambda2) - 2*alpha*(1/(lambda1*lambda1) + 1/(lambda2*lambda2))...
            -0.5 * epsilon * lambda1 * lambda1 * lambda2 * lambda2;
        if s2 > 0
            dydt = ode_fun1(t,y);
        else
            dydt = ode_fun2(t,y);
            
        end
    end

    function Mass = mass_matrix1(t,y)
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

    function Mass = mass_matrix2(t,y)
        R = 1 + gamma*cos(t);
        lambda_1 = 1/gamma * sqrt(y(2)^2 + y(4)^2);
        P_t = P*H/R_b;
        delta_sqrt = sqrt(P_t^2 * lambda_1^2 - 4*alpha*epsilon*lambda_1^4 ...
            + 16*alpha^2*lambda_1^4 + 32 * alpha*lambda_1^2 - 4*epsilon*lambda_1^2 +16);
        lambda_2_s_sq = (P_t * lambda_1 + delta_sqrt)/...
            (4*lambda_1 + 4*alpha*lambda_1^3 - epsilon * lambda_1^3);
        T1 = (P_t^2*lambda_1+32*alpha^2*lambda_1^3 - 8*epsilon*alpha*lambda_1^3+32*alpha*lambda_1-4*epsilon*lambda_1);
        A = T1/(delta_sqrt*(4*lambda_1+4*alpha*lambda_1^3 - epsilon*lambda_1^3));
        B = (P_t * lambda_1 + delta_sqrt)*(12*alpha*lambda_1^2 - 3*epsilon*lambda_1^2 + 4)/...
            (4*lambda_1+4*alpha*lambda_1^3 - epsilon*lambda_1^3)^2;
        
        dlambda_2_s_sq = P_t / (4*lambda_1+4*alpha*lambda_1^3 - epsilon*lambda_1^3)+A-B;
        
        dA = (P_t^2 + 96*alpha^2*lambda_1^2 - 24*epsilon*alpha*lambda_1^2+32*alpha-4*epsilon)/...
            (delta_sqrt*(4*lambda_1+4*alpha*lambda_1^3 - epsilon*lambda_1^3))-...
            T1^2/(delta_sqrt^3*(4*lambda_1+4*alpha*lambda_1^3 - epsilon*lambda_1^3))-...
            T1*(12*alpha*lambda_1^2-3*epsilon*lambda_1^2+4)/...
            (delta_sqrt*(4*lambda_1+4*alpha*lambda_1^3 - epsilon*lambda_1^3)^2);
        
        dB = (P_t*lambda_1 + delta_sqrt)*(lambda_1*(24*alpha-6*epsilon))/...
            (4*lambda_1+4*alpha*lambda_1^3 - epsilon*lambda_1^3)^2 + ...
            (P_t+T1/delta_sqrt)*(12*alpha*lambda_1^2 - 3*epsilon*lambda_1^2 + 4)/...
            (4*lambda_1+4*alpha*lambda_1^3 - epsilon*lambda_1^3)^2-...
            (P_t * lambda_1 + delta_sqrt)*(12*alpha*lambda_1^2 - 3*epsilon*lambda_1^2 + 4)^2/...
            (4*lambda_1+4*alpha*lambda_1^3 - epsilon*lambda_1^3)^3;
        
        ddlambda_2_s_sq = -P_t*(12*alpha*lambda_1^2-3*epsilon*lambda_1^2+4)/...
            (4*lambda_1 + 4*alpha*lambda_1^3-epsilon*lambda_1^3)^2 +dA-dB;

        Z = (2+6/(lambda_1^4 * lambda_2_s_sq)+2/(lambda_1^3 * lambda_2_s_sq^2)*dlambda_2_s_sq+...
            (2/(lambda_1^3 * lambda_2_s_sq^2) + 2/(lambda_1^2 * lambda_2_s_sq^3)*dlambda_2_s_sq)*dlambda_2_s_sq+...
            (1- 1/(lambda_1^2 * lambda_2_s_sq^2))*ddlambda_2_s_sq)+...
            alpha*(6/lambda_1^4 +2*lambda_1*dlambda_2_s_sq + 2*lambda_2_s_sq +...
            (2/dlambda_2_s_sq^3*dlambda_2_s_sq + 2*lambda_1)* dlambda_2_s_sq +...
            (lambda_1^2-1/dlambda_2_s_sq^2)*ddlambda_2_s_sq)+...
            epsilon/4*(2*lambda_2_s_sq + 4*lambda_1*dlambda_2_s_sq+...
            lambda_1^2*ddlambda_2_s_sq);

        Y = (2*lambda_1 + dlambda_2_s_sq - 2/(lambda_1^3 * lambda_2_s_sq) - 1/(lambda_1^2 * lambda_2_s_sq^2)*dlambda_2_s_sq)+...
            alpha*(-2/lambda_1^3 - dlambda_2_s_sq/lambda_2_s_sq^2 + lambda_1^2*dlambda_2_s_sq + 2*lambda_1*lambda_2_s_sq)+...
            epsilon/4*(2*lambda_1*lambda_2_s_sq + lambda_1^2*dlambda_2_s_sq);
        
        W_1 = y(2)/(gamma*sqrt(y(2)^2+y(4)^2));
        W_2 = y(4)/(gamma*sqrt(y(2)^2+y(4)^2));
        U_1 = 1/(gamma*sqrt(y(2)^2+y(4)^2)) - y(2)^2/(gamma*(y(2)^2+y(4)^2)^(3/2));
        U_2 = - (y(2)*y(4))/(gamma*(y(2)^2+y(4)^2)^(3/2));
        V_1 = - (y(2)*y(4))/(gamma*(y(2)^2+y(4)^2)^(3/2));
        V_2 = 1/(gamma*sqrt(y(2)^2+y(4)^2)) - y(4)^2/(gamma*(y(2)^2+y(4)^2)^(3/2));
        
        A22 = R*(U_1*Y+Z*W_1^2);
        A24 = R*(U_2*Y+Z*W_1*W_2);
        A42 = R*(V_1*Y+Z*W_1*W_2);
        A44 = R*(V_2*Y+Z*W_2^2);
        
        Mass = [1, 0, 0, 0;
            0, A22, 0, A24;
            0, 0, 1, 0;
            0, A42, 0, A44];
    end

% Define ODE:
    function dydt = ode_fun1(t,y)
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

    function dydt = ode_fun2(t,y)
        lambda_1 = 1/gamma * sqrt(y(2)^2 + y(4)^2);
        P_t = P*H/R_b;
        delta_sqrt = sqrt(P_t^2 * lambda_1^2 - 4*alpha*epsilon*lambda_1^4 ...
            + 16*alpha^2*lambda_1^4 + 32 * alpha*lambda_1^2 - 4*epsilon*lambda_1^2 +16);
        lambda_2_s_sq = (P_t * lambda_1 + delta_sqrt)/...
            (4*lambda_1 + 4*alpha*lambda_1^3 - epsilon * lambda_1^3);
        
        T1 = (P_t^2*lambda_1+32*alpha^2*lambda_1^3 - 8*epsilon*alpha*lambda_1^3+32*alpha*lambda_1-4*epsilon*lambda_1);
        A = T1/(delta_sqrt*(4*lambda_1+4*alpha*lambda_1^3 - epsilon*lambda_1^3));
        B = (P_t * lambda_1 + delta_sqrt)*(12*alpha*lambda_1^2 - 3*epsilon*lambda_1^2 + 4)/...
            (4*lambda_1+4*alpha*lambda_1^3 - epsilon*lambda_1^3)^2;
        
        dlambda_2_s_sq = -P_t * ( (8*alpha-2*epsilon)/((4*lambda_1+4*alpha*lambda_1^2 - epsilon*lambda_1^2)^2) -...
            (2*(8*alpha*lambda_1 - 2*epsilon*lambda_1)^2)/((4*lambda_1+4*alpha*lambda_1^2 - epsilon*lambda_1^2)^3)) +A-B;
        
        Y = (2*lambda_1 + dlambda_2_s_sq - 2/(lambda_1^3 * lambda_2_s_sq) - 1/(lambda_1^2 * lambda_2_s_sq^2)*dlambda_2_s_sq)+...
            alpha*(-2/lambda_1^3 - dlambda_2_s_sq/lambda_2_s_sq^2 + lambda_1^2*dlambda_2_s_sq + 2*lambda_1*lambda_2_s_sq)+...
            epsilon/4*(6*lambda_1^5*lambda_2_s_sq^3 + 3*lambda_1^6*lambda_2_s_sq^2*dlambda_2_s_sq);
        
        dydt_2 = sin(t)*y(2)/sqrt(y(2)^2+y(4)^2) * Y - P * y(1)*y(4)/gamma;
        dydt_4 = sin(t)*y(4)/sqrt(y(2)^2+y(4)^2) * Y + P * y(1)*y(2)/gamma;
        dydt = [y(2); dydt_2; y(4); dydt_4];
    end
end