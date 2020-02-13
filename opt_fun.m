function Sf = opt_fun(x)
    [drhoatpi, etaatpi, ~, ~] = IVP_solver(x);
    Sf = sqrt(drhoatpi^2 + etaatpi^2);
end

