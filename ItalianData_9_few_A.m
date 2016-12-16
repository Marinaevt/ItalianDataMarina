    a = importdata('1_125_R1.csv', ',');
    P1 = a(:,1);
    H_exp1 = a(:,2);
    t1 = a(:,3);
    s0 = a(1, 4);
    sf1 = a(1, 5);
    P1 = round(P1/0.25)*0.25;
    a = importdata('125_15_R1.csv', ',');
    P2 = a(:,1);
    H_exp2 = a(:,2);
    t2 = a(:,3);
    s0 = a(1, 4);
    sf2 = a(1, 5);
    P2 = round(P2/0.25)*0.25;
    a = importdata('125_15_R3.csv', ',');
    P3 = a(:,1);
    H_exp3 = a(:,2);
    t3 = a(:,3);
    s0 = a(1, 4);
    sf3 = a(1, 5);
    P3 = round(P3/0.25)*0.25;
    a = importdata('15_175_R2.csv', ',');
    P4 = a(:,1);
    H_exp4 = a(:,2);
    t4 = a(:,3);
    s0 = a(1, 4);
    sf4 = a(1, 5);
    P4 = round(P4/0.25)*0.25;
    A1 = 7.60110e+00;
    A2 = 8.57396e+01;
    A3 = -4.22663e+00;
    A4 = -3.29953e+01;
% A1 = 7.62738205788024;
% A2 =81.4475066902598;
% A3 = -4.54620442980213;
% A4 = -31.3610714688371;
    R0 = 22.5;
    Rho0 = 3;
%     Pstart = P(1);
%     Y = 1;
%     j = 1;
%     figure('Name', char(name))
%     H2 = [];
%     for i = 1:numel(P)
%         if P(i) ~= Pstart
%             [T, Y] = ode45(@(t, h)Testdhdt2_A1_A4(P(j), s0, Rho0, R0, h, sf, A1, A2, A3, A4), t(j:i), Y(end));
%             p45 = plot(T, Y, 'm', 'DisplayName', 'Starting values');
%             H2 = [H2; Y];
%             hold on
%             j = i;
%             Pstart = P(i);
%         end
%     end
     [X_min, Err] = fminsearchbnd(@(x)minsearcher7_few_A(H_exp1, P1, s0, Rho0, R0, t1, sf1, H_exp2, P2, t2, sf2, H_exp3, P3, t3, sf3, H_exp4, P4, t4, sf4, x(1), x(2), x(3), x(4)), [A1, A2, A3, A4], [0 0 inf inf], [inf inf 0 0], optimset('Display', 'iter', 'MaxFunEvals', 3000));   
     X_min
     Err
     Pstart = P1(1);
     Y(end) = 1;
     j = 1;
    H1 = [];
    for i = 1:numel(P1)
        if P1(i) ~= Pstart
            [T, Y] = ode45(@(t, h)Testdhdt2_A1_A4(P1(j), s0, Rho0, R0, h, sf1, X_min(1), X_min(2), X_min(3), X_min(4)), t1(j:i), Y(end));
            pkmn = plot(T, Y, 'r', 'DisplayName', 'Inverse analysis values');
            hold on
            H1 = [H1; Y];
            j = i;
            Pstart = P1(i);
        end
    end
%     mineu = minim2(H_exp, h, R0, Rho0, t, P)
    minode1 = minim2(H_exp1, H1, R0, Rho0, t1, P1)
     Pstart = P2(1);
     Y(end) = 1;
     j = 1;
    H1 = [];
    for i = 1:numel(P2)
        if P2(i) ~= Pstart
            [T, Y] = ode45(@(t, h)Testdhdt2_A1_A4(P2(j), s0, Rho0, R0, h, sf2, X_min(1), X_min(2), X_min(3), X_min(4)), t2(j:i), Y(end));
            pkmn = plot(T, Y, 'k', 'DisplayName', 'Inverse analysis values');
            hold on
            H1 = [H1; Y];
            j = i;
            Pstart = P2(i);
        end
    end
%     mineu = minim2(H_exp, h, R0, Rho0, t, P)
    minode2 = minim2(H_exp2, H1, R0, Rho0, t2, P2)
     Pstart = P3(1);
     Y(end) = 1;
     j = 1;
    H1 = [];
    for i = 1:numel(P3)
        if P3(i) ~= Pstart
            [T, Y] = ode45(@(t, h)Testdhdt2_A1_A4(P3(j), s0, Rho0, R0, h, sf3, X_min(1), X_min(2), X_min(3), X_min(4)), t3(j:i), Y(end));
            pkmn = plot(T, Y, 'c', 'DisplayName', 'Inverse analysis values');
            hold on
            H1 = [H1; Y];
            j = i;
            Pstart = P3(i);
        end
    end
%     mineu = minim2(H_exp, h, R0, Rho0, t, P)
    minode3 = minim2(H_exp3, H1, R0, Rho0, t3, P3)
    Pstart = P4(1);
    Y(end) = 1;
    j = 1;
    H1 = [];
    for i = 1:numel(P4)
        if P4(i) ~= Pstart
            [T, Y] = ode45(@(t, h)Testdhdt2_A1_A4(P4(j), s0, Rho0, R0, h, sf4, X_min(1), X_min(2), X_min(3), X_min(4)), t4(j:i), Y(end));
            pkmn = plot(T, Y, 'g', 'DisplayName', 'Inverse analysis values');
            hold on
            H1 = [H1; Y];
            j = i;
            Pstart = P4(i);
        end
    end
%     mineu = minim2(H_exp, h, R0, Rho0, t, P)
    minode4 = minim2(H_exp4, H1, R0, Rho0, t4, P4)
    hold on
    pexp = plot(t1, H_exp1, 'r--', 'DisplayName', 'Experimental data');
    pexp = plot(t2, H_exp2, 'k--', 'DisplayName', 'Experimental data');
    pexp = plot(t3, H_exp3, 'c--', 'DisplayName', 'Experimental data');
    pexp = plot(t4, H_exp4, 'g--', 'DisplayName', 'Experimental data');
    xlabel('Time')
    ylabel('Height')
%     legend([pexp p45 pkmn])
    axis([0 inf 0 35])
    hold off
