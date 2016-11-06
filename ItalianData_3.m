[filename, path, Filter] = uigetfile({'*.csv'}, 'Select the .CSV file' ,'MultiSelect', 'on');
if Filter == 0
    return
end
if size(char(filename(1)), 2) == 1
    filenum = 1;
else
    filenum = numel(filename);
end
for fileit = 1:filenum
    if filenum == 1
        name = fullfile(path, filename);
        a = importdata(name, ',');
    else
        name = fullfile(path, filename(fileit));
        a = importdata(char(name), ',');
    end
    a(isnan(a)) = 0;
    k = 1028.254506;
    m = 0.504606103;
    n = 0;
%     k = 7.007073092271015e+02;    
%     m = 0.456210908180900;    
%     n = 1.704630547644973e-04;
    B = 0.5 + 1/(1.54506*(m + 1)^3.2066);

    R0 = 22.5;
    Rho0 = 3;
    

%    t = a(:,1);
    P = a(:,1);
    P = round(P/0.25)*0.25;
    H_exp = a(:,2);
    t = a(:,3);
    s0 = a(:, 4);
    s0 = s0(1);
    num = size(a, 1);
    h = zeros(num, 1);
    h(1) = s0;

    dhdt = zeros(num, 1);
    x = zeros(num, 1);
    eps = zeros(num, 1);
    sigma = zeros(num, 1);
    dsdh_s = zeros(num, 1);
    r0 = zeros(num, 1);
    s_s0 = zeros(num, 1);
    SR = zeros(num, 1);
    Pstart = P(1);
    Y = 1;
    test = [];
    j = 1;
    figure('Name', char(name))
    %syms h_sym
    H2 = [];
    for i = 1:numel(P)
        if P(i) ~= Pstart
            %ttt = Testdhdt(k, n, m, P(j), s0, Rho0, R0, Y(end), h_sym);
            %tttt = matlabFunction(ttt);
            [T, Y] = ode45(@(t, h)Testdhdt(k, n, m, P(j), s0, Rho0, R0, h), t(j:i), Y(end));
            p45 = plot(T, Y, 'm', 'DisplayName', 'Runge-Kutta ode45');
            H2 = [H2; Y];
            hold on
            j = i;
            Pstart = P(i);
        end
    end
    %p45 = plot(t(j:i), Y, 'm', 'DisplayName', 'Runge-Kutta ode45');
%     ttt = Testdhdt(k, n, m, P(1), s0, B, Rho0, R0, 1);
%     tttt = matlabFunction(ttt);
%     [T, Y] = ode45(@(t, h)tttt(h), t(1:433), 1);
%     plot(t(1:433), Y, 'k')
%     hold on
%     ttt = Testdhdt(k, n, m, 1.5, s0, B, Rho0, R0, 1);
%     tttt = matlabFunction(ttt);
%     [T, Y1] = ode45(@(t, h)tttt(h), t(434:583), Y(end));
%     plot(t(434:583), Y1, 'k')
%     ttt = Testdhdt(k, n, m, 1, s0, B, Rho0, R0, 1);
%     tttt = matlabFunction(ttt);
%     [T, Y2] = ode45(@(t, h)tttt(h), t(584:735), Y1(end));
%     plot(t(584:735), Y2, 'k')

    for i = 1:size(P)
        if i > 1
            h(i) = h(i - 1) + dhdt(i) * (t(i) - t(i - 1));
        else
            h(i) = 1;
        end
        if h(i) < R0 + Rho0
            r0(i) = ((R0 + Rho0)^2 + h(i)^2)/2/h(i) - Rho0;
        else
            r0(i) = R0;
        end
        x(i) = min(h(i)/(r0(i) + Rho0), 1);
        s_s0(i) = 1 - B * x(i);
        sigma(i) = P(i) * r0(i) / (s_s0(i) * s0 * 2);
        eps(i) = -log(s_s0(i));
        dsdh_s(i) = -B*(h(i) - 2*(r0(i) + Rho0))/((r0(i) + Rho0)^2) / s_s0(i);
        if x(i) < 0.999% <= 0.01
            dhdt(i + 1) = dHdT(sigma(i), k, eps(i), n, m, dsdh_s(i));
        else
            dhdt(i + 1) = 0.1;
        end
        if i > 1
            SR(i) = (eps(i) - eps(i-1)) / (t(i) - t(i-1));
        end
    end
    SR(1) = SR(2);
    F = [];
    F1 = [];
%     num0 = find(H_exp == 0);
%     num0 = num0(1) - 1;
%     num10 = find(t == 10)-1;
%     pit = P(num10+1);
    k = 1;
%    syms H_exp_sym k_sym n_sym m_sym
%     test = Testdhdt(k_sym, n_sym, m_sym, P(1+num10), s0, B, Rho0, R0, h(1+num10), h_sym);
%     testt = minim1(H_exp_sym, test);
%     testfun = matlabFunction(testt);
%     [X_min, Err_min] = fminsearch(@(x)testfun(H_exp(1), h(1+num10), x(3), x(4), x(5)), [k, m, n]);

     [X_min, Err] = fminsearchbnd(@(x)minsearcher1(H_exp, x(1), x(2), x(3), P, s0, Rho0, R0, t), [k, m, n], [100 inf inf], [10000 inf inf], optimset('Display', 'iter', 'MaxFunEvals', 3000));
%    [X_min, Err] = fminsearch(@(x)minsearcher1(H_exp, x(1), x(2), x(3), P, s0, Rho0, R0, t, num0), [k, m, n], optimset('Display', 'iter'));
    
     X_min
     Err
     Pstart = P(1);
    Y(end) = 1;
     j = 1;
%     hold on
%     %figure('Name', char(name))
%     %syms h_sym
    H1 = [];
    for i = 1:numel(P)
        if P(i) ~= Pstart
            %ttt = Testdhdt(k, n, m, P(j), s0, Rho0, R0, Y(end), h_sym);
            %tttt = matlabFunction(ttt);
            [T, Y] = ode45(@(t, h)Testdhdt(X_min(1), X_min(3), X_min(2), P(j), s0, Rho0, R0, h), t(j:i), Y(end));
            pkmn = plot(T, Y, 'r', 'DisplayName', 'k, m and n from fminsearch');
            hold on
            H1 = [H1; Y];
            j = i;
            Pstart = P(i);
        end
    end
    mineu = minim2(H_exp, h, R0, Rho0, t, P)
    minode = minim2(H_exp, H1, R0, Rho0, t, P)


    pcalc = plot(t, h, 'b--', 'DisplayName', 'Euler method');
    hold on
%      plot(F(:, 3), F(:, 1), 'r', 'DisplayName', 'Data considering error')
    pexp = plot(t, H_exp, 'g', 'DisplayName', 'Experimental data');
    xlabel('Time')
    ylabel('Height')
%     legend([pcalc pexp p45 pkmn])
    axis([0 inf 0 35])
    hold off
%     figure('Name', char(name))
%     err = MSE(H_exp(1:num0), h(num10:numel(h)));
%     err(1:numel(t)) = err;
%     pmse = plot(t, err, 'r', 'DisplayName', 'MSE');
%     hold on
%     err(1:numel(t)) = err.^2;
%     pdx = plot(t, err, 'y', 'DisplayName', 'DX');
%     xlabel('Time')
%     ylabel('Error')
%     legend([pmse pdx])
% %     plot(t(1:num0), F1(:, 1), 'y', 'DisplayName', 'Data considering error ^2')
% %    plot(t, H_exp, 'g', 'DisplayName', 'Experimental data');
%     hold off
%     figure('Name', char(name))
%     plot(t, SR)
%     xlabel('Time')
%     ylabel('SR')
%     axis([0 inf 0 0.01])
end
