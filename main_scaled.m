clc,clear all;

% % Global avariables:
global epsilon alpha gamma t_span H rho0 scale_coeff R_b

%% thickness
H = 0.01;
%% radius
R_b = 1;
%% electrical load : epsilon = phi_0^2 / (C_1 * beta * H^2);
epsilon = 0.;
%% alpha = C_2/C_1 (material property)
alpha = 0.3;
%% aspect ratio of radii (geometry parameter)
gamma = 0.6;
%% gamma = 0.4 alpha = 0.2 E = 0
% there are some more initial guesses for different gamma and alpha values in the
% 'initial_guess.txt' file.
% rho0 = 1.49;
% detai0 = 0.4691;
% P0 = 3.5424;
%% gamma = 0.6 alpha = 0.3 E = 0.
% rho0 = 1.7;
% detai0 = 0.6813;
% P0 = 2.2933;
rho0 = 3.225;
detai0 = 1.4118;
P0 = 4.005;
%%
warning off;
% step size of rho
rho_step = 0.01;
% number of steps
n_step = 100;
% range of theta:
n_elements = 360;
t_span = linspace(0,pi,n_elements+1);
 % initialise some vectors
rho0_vec = zeros(n_step,1);
v_ratio = zeros(n_step,1);
Pvec=zeros(n_step,1);
%% An initial guess of scaling coefficient to make rho_step/P_step close to 1
%scale_coeff = 1;
scale_coeff = 0.05;
% rho0 = rho0 + rho_step + rho_step2;
for i = 1:n_step
% rho0 = rho0-rho_step2
rho0 = rho0+rho_step;
 if i > 2
    P_step = Pvec(i-2)-Pvec(i-1);
    scale_coeff = abs(rho_step/P_step);
 end
 rho0_vec(i,1) = rho0;
% scale x0 using the coefficient
 x0_s=[P0,detai0/scale_coeff];

 options = optimset('FunValCheck', 'on', 'MaxFunEvals', 10000, 'TolFun', 1e-14, 'TolX', 1e-6);
 
 [x, fval, exitflag, output] = fminsearchbnd(@opt_fun_2, x0_s, [0,0], [ ], options);
 
 [drhoatpi, etaatpi, t, y, dy] = IVP_solver_with_derivatives(x);
 
 [s11,s22,is_negative,v_negative_stress,p_negative_stress] = principle_stresses(y,x,t);
 
 figure(1)
 grid on;
 plot(y(:,1), y(:,3), 'linewidth', 0.5);
 hold on;
 
 % check second variation for wavenumber = 1,2,3,4
 for wavenumber = 1:4
     [l1,l2,detM,iszero,P_bifurcation,v_ratio_bifurcation] = check_second_variation( x, y, dy, wavenumber);
     if iszero == 1
         disp("With wave number = " + string(wavenumber) + ", we have bifurcation when volume change = "...
             +string(v_ratio_bifurcation(1)) + " and critial load = "+ string(P_bifurcation(1)));
     end
     figure(2)
     grid on;
     plot(t_span,detM,'-b')
     hold on;
 end
 
 if is_negative == 1
     % here we use the results as the initial guess for solving the problem
     % with energy relaxation in the region will S_22 is negative.
     P0 = x(1);
     detai0 = y(1,4);
     x0=[P0,detai0];
     [x, fval, exitflag, output] = fminsearchbnd(@opt_fun_ER, x0, [0,0], [], options);
     [drhoatpi2, etaatpi2, t2, y] = IVP_solver_ER(x);
 end
 
 v_ratio(i,1) = volume_change(y);
 
 Pvec(i,1) = x(1);
 
 P0 = x(1);
 detai0 = y(1,4);
end
 figure(5)
 plot(v_ratio,Pvec, 'linewidth', 1.5);
 hold on;
 a = [v_ratio'; Pvec'];
 
 %% gamma = 0.4 alpha = 0.2 E = 0
rho0 = 1.49 + rho_step;
detai0 = 0.4691;
P0 = 3.5424;

% step size of rho
rho_step_backward = 0.0025;
% number of steps
n_step_backward = 20;
 % initialise some vectors
rho0_vec = zeros(n_step_backward,1);
v_ratio = zeros(n_step_backward,1);
Pvec=zeros(n_step_backward,1);
%%
for i = 1:n_step_backward
    
rho0 = rho0-rho_step_backward

 if i > 2
    P_step = Pvec(i-2)-Pvec(i-1);
    scale_coeff = abs(n_step_backward/P_step);
 end

 rho0_vec(i,1) = rho0;
 
 x0_s=[P0,detai0/scale_coeff];

 options = optimset('FunValCheck', 'on', 'MaxFunEvals', 10000, 'TolFun', 1e-14, 'TolX', 1e-6);
 
 [x, fval, exitflag, output] = fminsearchbnd(@opt_fun_2, x0_s, [0,0], [ ], options);
 
 [drhoatpi, etaatpi, t, y, dy] = IVP_solver_with_derivatives(x);
 
 [s11,s22,is_negative,v_negative_stress,p_negative_stress] = principle_stresses(y,x,t);
 
 figure(1)
 grid on;
 plot(y(:,1), y(:,3), 'linewidth', 0.5);
 hold on;
 
 v_ratio(i,1) = volume_change(y);
 
 Pvec(i,1) = x(1);
 
 P0 = x(1);
 detai0 = y(1,4);
end
 figure(5)
 plot(v_ratio,Pvec, 'linewidth', 1.5);
 hold on;
 b = [v_ratio'; Pvec'];
 b = sort(b,2);
 a = [b a];
 fileID = fopen("gamma = "+string(gamma)+ ", alpha = "+string(alpha) + ...
     ", epsilon = "+ string(epsilon) + ' .txt','w');
 formatSpec = '%4.6f %4.3f\n';
 fprintf(fileID,formatSpec,a);
 fclose(fileID);