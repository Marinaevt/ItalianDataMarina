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
    B = 0.5 + 1/(1.54506*(m + 1)^3.2066);

    R0 = 22.5;
    Rho0 = 3;
    s0 = 1;

    t = a(:,1);
    P = a(:,2);
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
    j = 1;
    figure('Name', char(name))
    %syms h_sym
    for i = 1:numel(P)
        if P(i) ~= Pstart
            %ttt = Testdhdt(k, n, m, P(j), s0, Rho0, R0, Y(end), h_sym);
            %tttt = matlabFunction(ttt);
            [T, Y] = ode45(@(t, h)Testdhdt(k, n, m, P(j), s0, Rho0, R0, h), t(j:i), Y(end));
            p45 = plot(t(j:i), Y, 'm', 'DisplayName', 'Runge-Kutta ode45');
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
%             ttt = Testdhdt(k, n, m, P(i), s0, B, Rho0, R0, h(i-1));
%             tttt = matlabFunction(ttt);
%             [T, Y] = ode45(@(t, h)tttt(h), [t(1), t(end)], h(i-1));
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
    H_exp = a(:,3);
    t_exp = a(:,4);
    F = [];
    F1 = [];
    num0 = find(H_exp == 0);
    num0 = num0(1) - 1;
    num10 = find(t == 10)-1;
    pit = P(num10+1);
    k = 1;
    syms H_exp_sym k_sym n_sym m_sym
%     test = Testdhdt(k_sym, n_sym, m_sym, P(1+num10), s0, B, Rho0, R0, h(1+num10), h_sym);
%     testt = minim1(H_exp_sym, test);
%     testfun = matlabFunction(testt);
%     [X_min, Err_min] = fminsearch(@(x)testfun(H_exp(1), h(1+num10), x(3), x(4), x(5)), [k, m, n]);

    

    Err_min = 10000;
    Pstart = P(1+num10);
    %wait = waitbar(0, '1', 'Name', 'Calculating Error');
    Y_end = H_exp(1);
    j = 1;
    flag = 0;
    kmn = [];
    test = [];
    for i = 1:num0
        if P(i+num10) ~= Pstart
            %ttt = Testdhdt(k, n, m, P(j+num10), s0, Rho0, R0, Y_end, h_sym);
            %tttt = matlabFunction(ttt);
            if flag == 0
                [T, Y] = ode45(@(t, h)Testdhdt(k, n, m, P(j), s0, Rho0, R0, h), t(j+num10:i+num10), Y_end);
                Err_min = minim3(H_exp(j:i), Y, R0);
                flag = 1;
            end
%             x(1) = k
%             x(2) = m
%             x(3) = n
%            test = Testdhdt(k_sym, n_sym, m_sym, P(i+num10), s0, B, Rho0, R0, h(i+num10), h_sym);
%            testfun = matlabFunc(testt);
            [X, Err] = fminsearchbnd(@(x)minsearcher(H_exp(j:i), x(1), x(2), x(3), P(j), s0, Rho0, R0, Y_end, t(j+num10:i+num10)), [k, m, n], [100 inf inf], [30000 inf inf], optimset('Display', 'iter'));
            if Err_min >= Err
                Err_min = Err;
                X_min = X;
            end
            kmn = [kmn; X];
            test = [test; minsearcher1(H_exp, X(1), X(2), X(3), P, s0, Rho0, R0, t, num0)];
            Y_end = Y(end);
            Pstart = P(i+num10);
            j = i;
        end
        %waitbar(i / num0, wait, sprintf('%d/%d points calculated', i, num0), 'r');
    end
    %close(wait)
    Err_min
    X_min
    test
    kmn
    Pstart = P(1);
    Y(end) = 1;
    j = 1;
    hold on
    %figure('Name', char(name))
    %syms h_sym
    for i = 1:numel(P)
        if P(i) ~= Pstart
            %ttt = Testdhdt(k, n, m, P(j), s0, Rho0, R0, Y(end), h_sym);
            %tttt = matlabFunction(ttt);
            [T, Y] = ode45(@(t, h)Testdhdt(X_min(1), X_min(3), X_min(2), P(j), s0, Rho0, R0, h), t(j:i), Y(end));
            p45 = plot(t(j:i), Y, 'r', 'DisplayName', 'Runge-Kutta ode45');
            hold on
            j = i;
            Pstart = P(i);
        end
    end
%    wisdom1 = @(x)(((x(2) - x(1))/x(2)).^2 + ((x(4) - x(3))/x(4)).^2); 
%     wait = waitbar(0, '1', 'Name', 'Calculating Error');
%     for i = 1:num0 % 3 столбец - время, 4 - время (эксп.)
%         Fmin = fminsearch(wisdom1,[h(i+num10), H_exp(i), t(i+num10), t_exp(i)]);
%         Fmin1 = fminsearch(@wisdom,[h(i+num10), H_exp(i)]);
%         F = [F; Fmin];
%         F1 = [F1; Fmin1];
%         waitbar(i / num0, wait, sprintf('%d/%d points calculated', i, num0), 'r')
%     end
%     close(wait)
%     figure('Name', char(name))
    pcalc = plot(t, h, 'b--', 'DisplayName', 'Euler method');
    %hold on
%      plot(F(:, 3), F(:, 1), 'r', 'DisplayName', 'Data considering error')
    pexp = plot(t_exp(1:num0), H_exp(1:num0), 'g', 'DisplayName', 'Experimental data');
    xlabel('Time')
    ylabel('Height')
    legend([pcalc pexp p45])
    axis([0 inf 0 35])
    hold off
    figure('Name', char(name))
    err = MSE(H_exp(1:num0), h(num10:numel(h)));
    err(1:numel(t)) = err;
    pmse = plot(t, err, 'r', 'DisplayName', 'MSE');
    hold on
    err(1:numel(t)) = err.^2;
    pdx = plot(t, err, 'y', 'DisplayName', 'DX');
    xlabel('Time')
    ylabel('Error')
    legend([pmse pdx])
%     plot(t_exp(1:num0), F1(:, 1), 'y', 'DisplayName', 'Data considering error ^2')
%    plot(t_exp, H_exp, 'g', 'DisplayName', 'Experimental data');
    hold off
    figure('Name', char(name))
    plot(t, SR)
    xlabel('Time')
    ylabel('SR')
    axis([0 inf 0 0.01])
end
