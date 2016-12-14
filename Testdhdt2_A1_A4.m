function [dhdt] = Testdhdt2(P, s0, Rho0, R0, h, sf, A1, A2, A3, A4)
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
    
    dsdh = B * (h - 2*(r0 + Rho0))/((r0 + Rho0)^2) * s0;
%     dsdh = (4*B*h^3*s0)/((R0 + Rho0)^2 + h^2)^2 - (4*B*h*s0)/((R0 + Rho0)^2 + h^2);
%     dsdh = (2*s0*exp(-(4*Rho0*(atan(Rho0/(R0^2 + 2*Rho0*R0)^(1/2)) - atan((Rho0 - h)/(R0^2 + 2*Rho0*R0)^(1/2))))/(R0^2 + 2*Rho0*R0)^(1/2))*(R0 + Rho0)^4*(2*Rho0 - 2*h))/((R0 + Rho0)^2 - 2*Rho0*h + h^2)^3 - (4*Rho0*s0*exp(-(4*Rho0*(atan(Rho0/(R0^2 + 2*Rho0*R0)^(1/2)) - atan((Rho0 - h)/(R0^2 + 2*Rho0*R0)^(1/2))))/(R0^2 + 2*Rho0*R0)^(1/2))*(R0 + Rho0)^4)/(((Rho0 - h)^2/(R0^2 + 2*Rho0*R0) + 1)*(R0^2 + 2*Rho0*R0)*((R0 + Rho0)^2 - 2*Rho0*h + h^2)^2);
%     dsdh = -2*s/r0;
%     k = 583.9272;
%     m = 0.5418;
%     k = 1028.254506;
%     m = 0.504606103;

% sref = 0.0000220769518118314;
% s_0 = 4.28166914057785;
% s_1 = 12.98367869;
% m0 = 0.654006407143157;
% mju = -0.54197793;
% m = m0 + mju*eps;
% k = (s_0 + s_1*eps)/sref^(m);
    if x < 0.999
%         stress = A1*ln(strain_rate)+A2 + strain*(A3*ln(strain_rate)+A4);
        eps_dot = exp((sigma - eps*A4 - A2)/(A1+eps*A3));
%         eps_dot = (sigma/k)^(1/m);
        dhdt = max(0, min(1000000, -s/dsdh*eps_dot));
    else
        dhdt = 0.1;
    end
end

