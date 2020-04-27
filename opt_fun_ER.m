function Sf = opt_fun_ER(x)
% cost function with energy relaxation
    [drhoatpi, etaatpi, ~, ~] = IVP_solver_ER(x);
    Sf = sqrt(drhoatpi^2 + etaatpi^2);
end

