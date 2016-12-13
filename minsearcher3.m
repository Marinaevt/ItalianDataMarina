function [fun] = minsearcher3(H_exp, P, s0, Rho0, R0, t, sf, A1, A2, A3, A4)
    Pstart = P(1);
    Y = 1;
    j = 1;
    H = [];
    for i = 1:numel(P)
        if P(i) ~= Pstart
            [T, Y] = ode45(@(t, h)Testdhdt2(P(j), s0, Rho0, R0, h, sf, A1, A2, A3, A4), t(j:i), Y(end));
            H = [H; Y];
            j = i;
            Pstart = P(i);
        end
    end
    fun = minim2(H_exp, H, R0, Rho0, t, P);
end

