function [fun] = minsearcher(H_exp, k, m, n, P, s0, Rho0, R0, Y, t)
    %ttt = Testdhdt(k, n, m, P, s0, Rho0, R0, Y, h_sym);
    %tttt = matlabFunction(ttt);
    [T, H] = ode45(@(t, h)Testdhdt(k, n, m, P, s0, Rho0, R0, h), t, Y);
    fun = minim3(H_exp, H, R0);
end

