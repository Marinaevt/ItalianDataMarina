function [dhdt] = Testdhdt3(k, n, m, P, s0, Rho0, R0, h)
    B = 0.5 + 1/(1.54506*(m + 1)^3.2066);
    if h < R0 + Rho0
        r0 = ((R0 + Rho0)^2 + h^2)/2/h - Rho0;
    else
        r0 = R0;
    end
    x = min(h/(r0 + Rho0), 1);
    s_s0 = 1 - B * x;
    sigma = P * r0 / (s_s0 * s0 * 2);
    eps = -log(s_s0);
    dsdh_s = -B*(h - 2*(r0 + Rho0))/((r0 + Rho0)^2) / s_s0;
    if x < 0.999
        dhdt = max(0, min(1000000, (sigma / (k * eps^n))^(1/m))/(dsdh_s));
    else
        dhdt = 0.1;
    end
    %dhdt = (sigma / (k * eps^n))^(1/m)/(dsdh_s);
end

