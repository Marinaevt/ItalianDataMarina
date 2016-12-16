function [dhdt, eps_dot] = Testdhdt4_km_s0s1(P, s0, Rho0, R0, h, sf, t, eps_dot, s_0, s_1)
    if h < R0 + Rho0
        r0 = ((R0 + Rho0)^2 + h^2)/2/h - Rho0;
    else
        r0 = R0;
    end
    x = min(h/(r0 + Rho0), 1);
    B = 1 - sf/s0;  
    s = s0 * (1 - B*(h/(r0 + Rho0)));
    sigma = P * r0 / (s * 2);
    eps = log(s0/s);
    A1 = 7.60110e+00;
    A2 = 8.57396e+01;
    A3 = -4.22663e+00;
    A4 = -3.29953e+01;
    dsdh = B * (h - 2*(r0 + Rho0))/((r0 + Rho0)^2) * s0;
    m = 0.627 - t*1.02367e-5;
%     for i = 1:120
%         sigma1 =  A1*log(eps_dot)+A2 + eps*(A3*log(eps_dot)+A4);
%         k = sigma1/eps_dot^m;
%         eps_dot = (sigma/k)^(1/m);
%     end
	  sigma1 =  s_0 + eps*s_1;
%     sigma1 =  4.2817 + eps*12.9837;
%     sigma1 =  4.2817 + t*0.000286640048804879;
    k = sigma1/(0.000340154387908685)^m;
    eps_dot = (sigma/k)^(1/m);

    if x < 0.999
%         stress = A1*ln(strain_rate)+A2 + strain*(A3*ln(strain_rate)+A4);
%         eps_dot = exp((sigma - eps*A4 - A2)/(A1+eps*A3));
        eps_dot = (sigma/k)^(1/m);
        dhdt = max(0, min(1000000, -s/dsdh*eps_dot));
    else
        dhdt = 0.1;
    end
end

