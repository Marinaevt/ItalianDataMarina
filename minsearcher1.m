function [fun] = minsearcher1(H_exp, k, m, n, P, s0, Rho0, R0, t, num0)
    Pstart = P(1);
    Y = 1;
    j = 1;
    H = [];
    for i = 1:numel(P)
        if P(i) ~= Pstart
            [T, Y] = ode45(@(t, h)Testdhdt(k, n, m, P(j), s0, Rho0, R0, h), t(j:i), Y(end));
            H = [H; Y];
            j = i;
%             plot(t(j:i), Y);
%             hold on
            Pstart = P(i);
        end
    end
%     k
%     m
%     n
    fun = minim2(H_exp, H, num0, R0, Rho0, t, P);
end

