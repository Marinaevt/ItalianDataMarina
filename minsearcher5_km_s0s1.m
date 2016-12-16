function [fun] = minsearcher5_km_s0s1(H_exp, P, s0, Rho0, R0, t, sf, s_0, s_1)
    Pstart = P(1);
    Y = 1;
    j = 1;
    H = [];
    eps_dot = 0.0001;
    for i = 1:numel(P)
        if P(i) ~= Pstart
            [T, Y] = ode45(@(t, h)Testdhdt4_km_s0s1(P(j), s0, Rho0, R0, h, sf, t, eps_dot, s_0, s_1), t(j:i), Y(end));
            H = [H; Y];
            j = i;
            Pstart = P(i);
        end
    end
    fun = minim2(H_exp, H, R0, Rho0, t, P);
end

