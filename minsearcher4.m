function [fun] = minsearcher4(H_exp, P, s0, Rho0, R0, t, k0, m0, n, mju, psy)
    Pstart = P(1);
    Y = 1;
    j = 1;
    H = [];
    for i = 1:numel(P)
        if P(i) ~= Pstart 
%             ttest = t(j:i);
%             test = t(i);
            m = m0 + mju*t(j);
            k = k0 + psy*t(j);
            [T, Y] = ode45(@(t, h)Testdhdt3(k, n, m, P(j), s0, Rho0, R0, h), t(j:i), Y(end));
            H = [H; Y];
            j = i;
            Pstart = P(i);
%             k0 = k;
%             m0 = m;
        end
    end
    fun = minim2(H_exp, H, R0, Rho0, t, P);
end

