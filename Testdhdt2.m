function [dhdt] = Testdhdt2(P, s0, Rho0, R0, h, A1, A2, A3, A4)
%     B = 0.5 + 1/(1.54506*(m + 1)^3.2066);
    if h < R0 + Rho0
        r0 = ((R0 + Rho0)^2 + h^2)/2/h - Rho0;
    else
        r0 = R0;
    end
     x = min(h/(r0 + Rho0), 1);
%     s_s0 = 1 - B * x; %??????????
    C = (R0^2 + 2*R0*Rho0)^(-1/2);
    s = (s0*(R0+Rho0)^4)/((R0+Rho0)^2 + h^2 - 2*h*Rho0)^2 * exp(-4*C*Rho0*(atan(C*Rho0) - atan(C*(Rho0 - h))));
    sigma = P * r0 / (s * 2);
    eps = -log(s0/s);
%     dsdh_s = -B*(h - 2*(r0 + Rho0))/((r0 + Rho0)^2) / s_s0;
    if x < 0.999
        %stress = A1*ln(strain_rate)+A2 + strain*(A3*ln(strain_rate)+A4)
        dhdt = max(0, min(1000000, 1/2 * exp((sigma - eps*A4 - A2)/(A1+eps*A3))*r0));
    else
        dhdt = 0.1;
    end
    %dhdt = (sigma / (k * eps^n))^(1/m)/(dsdh_s);
end

