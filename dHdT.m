function [dhdt] = dHdT(sigma, k, eps, n, m, dsdh_s)
    dhdt = max(0, min(1000000, (sigma / (k * eps^n))^(1/m))/(dsdh_s));
end

