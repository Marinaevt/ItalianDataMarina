function [fun] = minsearcher6(H_exp, P, s0, Rho0, R0, t, sf, sref, s_0, s_1, m0, mju)
    Pstart = P(1);
    Y = 1;
    j = 1;
    H = [];
    for i = 1:numel(P)
        if P(i) ~= Pstart
            [T, Y] = ode45(@(t, h)Testdhdt5(P(j), s0, Rho0, R0, h, sf, sref, s_0, s_1, m0, mju), t(j:i), Y(end));
            H = [H; Y];
            j = i;
            Pstart = P(i);
        end
    end
    fun = minim2(H_exp, H, R0, Rho0, t, P);
end

