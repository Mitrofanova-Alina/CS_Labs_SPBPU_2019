function [result] = IVE(A_inf, A_sup, b_inf, b_sup, tolmax, arg, n)

    mid_b = (b_inf + b_sup) / 2;
    rad_b = (b_inf - b_sup) / 2;
    b_appox = 1/2 * (abs(mid_b +rad_b) + abs(mid_b - rad_b));

    minCond = HeurMinCond(A_inf, A_sup);
    arg_norm = norm(arg);
    b_norm = norm(b_appox);

    result = sqrt(n) * tolmax * minCond * arg_norm  / b_norm ;

end