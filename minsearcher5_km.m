function [fun] = minsearcher5_km(H_exp, P, s0, Rho0, R0, t, sf, k, m)
    Pstart = P(1);
    Y = 1;
    j = 1;
    H = [];
    for i = 1:numel(P)
        if P(i) ~= Pstart
            [T, Y] = ode45(@(t, h)Testdhdt4_km(P(j), s0, Rho0, R0, h, sf, k, m), t(j:i), Y(end));
            H = [H; Y];
            j = i;
            Pstart = P(i);
        end
    end
    fun = minim2(H_exp, H, R0, Rho0, t, P);
end

