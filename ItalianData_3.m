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

     [X_min, Err] = fminsearchbnd(@(x)minsearcher1(H_exp, x(1), x(2), x(3), P, s0, Rho0, R0, t), [k, m, n], [100 inf inf], [10000 inf inf], optimset('Display', 'iter', 'MaxFunEvals', 3000));
    
     X_min
     Err
     Pstart = P(1);
    Y(end) = 1;
     j = 1;
    H1 = [];
    for i = 1:numel(P)
        if P(i) ~= Pstart
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
    pexp = plot(t, H_exp, 'g', 'DisplayName', 'Experimental data');
    xlabel('Time')
    ylabel('Height')
%     legend([pcalc pexp p45 pkmn])
    axis([0 inf 0 35])
    hold off
end
