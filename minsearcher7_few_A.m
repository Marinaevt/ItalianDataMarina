function [fun] = minsearcher7(H_exp1, P1, s0, Rho0, R0, t1, sf1, H_exp2, P2, t2, sf2, H_exp3, P3, t3, sf3, H_exp4, P4, t4, sf4, A1, A2, A3, A4)
H1 = [];        
Pstart = P1(1);
        Y = 1;
        j = 1;
        H = [];
        for i = 1:numel(P1)
            if P1(i) ~= Pstart
                [T, Y] = ode45(@(t, h)Testdhdt2_A1_A4(P1(j), s0, Rho0, R0, h, sf1, A1, A2, A3, A4), t1(j:i), Y(end));
                H = [H; Y];
                j = i;
                Pstart = P1(i);
            end
        end
        H1 = [H1; H];
    %fun1 = minim2(H_exp1, H, R0, Rho0, t1, P1);
        Pstart = P2(1);
        Y = 1;
        j = 1;
        H = [];
        for i = 1:numel(P2)
            if P2(i) ~= Pstart
                [T, Y] = ode45(@(t, h)Testdhdt2_A1_A4(P2(j), s0, Rho0, R0, h, sf2, A1, A2, A3, A4), t2(j:i), Y(end));
                H = [H; Y];
                j = i;
                Pstart = P2(i);
            end
        end
          H1 = [H1; H];
    %fun2 = minim2(H_exp2, H, R0, Rho0, t2, P2);
        Pstart = P3(1);
        Y = 1;
        j = 1;
        H = [];
        for i = 1:numel(P3)
            if P3(i) ~= Pstart
                [T, Y] = ode45(@(t, h)Testdhdt2_A1_A4(P3(j), s0, Rho0, R0, h, sf3, A1, A2, A3, A4), t3(j:i), Y(end));
                H = [H; Y];
                j = i;
                Pstart = P3(i);
            end
        end
          H1 = [H1; H];
        %     fun3 = minim2(H_exp3, H, R0, Rho0, t3, P3);
            Pstart = P4(1);
        Y = 1;
        j = 1;
        H = [];
        for i = 1:numel(P4)
            if P4(i) ~= Pstart
                [T, Y] = ode45(@(t, h)Testdhdt2_A1_A4(P4(j), s0, Rho0, R0, h, sf4, A1, A2, A3, A4), t4(j:i), Y(end));
                H = [H; Y];
                j = i;
                Pstart = P4(i);
            end
        end
          H1 = [H1; H];
          H_exp_all = [H_exp1; H_exp2; H_exp3; H_exp4];
%     fun4 = minim2(H_exp4, H, R0, Rho0, t4, P4);
%    fun = (fun1+fun2+fun3+fun4)/4;
fun = minim2(H_exp_all, H1, R0, Rho0, t1, P4);
end

