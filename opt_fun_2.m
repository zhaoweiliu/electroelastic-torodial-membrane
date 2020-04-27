function Sf = opt_fun_2(x)
% cost function with scaled values
    [drhoatpi, etaatpi, ~, ~] = IVP_solver_scaled(x);
    Sf = sqrt(drhoatpi^2 + etaatpi^2);
end