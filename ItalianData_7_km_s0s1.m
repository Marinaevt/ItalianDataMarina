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
    k = 1028.254506;
    m = 0.504606103;
    R0 = 22.5;
    Rho0 = 3;
    
    P = a(:,1);
    P = round(P/0.25)*0.25;
    H_exp = a(:,2);
     H_exp = find(H_exp>5);
    t = a(:,3);
    s0 = a(:, 4);
    s0 = s0(1);
    sf = a(:, 5);
    sf = sf(1);
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
    H2 = [];
    eps_dot = 0.00001;
    s_0 = 25.0337982809274;
    s_1 = 0.678009796337607;
    for i = 1:numel(P)
        if P(i) ~= Pstart
            [T, Y] = ode45(@(t, h)Testdhdt4_km_s0s1(P(j), s0, Rho0, R0, h, sf, t, eps_dot, s_0, s_1), t(j:i), Y(end));
            p45 = plot(T, Y, 'm', 'DisplayName', 'Starting values');
            H2 = [H2; Y];
            hold on
            j = i;
            Pstart = P(i);
        end
    end

     [X_min, Err] = fminsearchbnd(@(x)minsearcher5_km_s0s1(H_exp, P, s0, Rho0, R0, t, sf, x(1), x(2)), [s_0, s_1], [inf inf], [inf inf], optimset('Display', 'iter', 'MaxFunEvals', 3000));   
     X_min
     Err
     Pstart = P(1);
     Y(end) = 1;
     j = 1;
    H1 = [];
    for i = 1:numel(P)
        if P(i) ~= Pstart
            [T, Y] = ode45(@(t, h)Testdhdt4_km_s0s1(P(j), s0, Rho0, R0, h, sf, t, eps_dot, X_min(1), X_min(2)), t(j:i), Y(end));
            pkmn = plot(T, Y, 'r', 'DisplayName', 'Inverse analysis values');
            hold on
            H1 = [H1; Y];
            j = i;
            Pstart = P(i);
        end
    end
    minode = minim2(H_exp, H1, R0, Rho0, t, P)

    hold on
    pexp = plot(t, H_exp, 'g', 'DisplayName', 'Experimental data');
    xlabel('Time')
    ylabel('Height')
    legend([pexp p45 pkmn])
    axis([0 inf 0 35])
    savefig(strcat(strrep(filename, '.csv', ''), '_class_Bakofen.fig'))
    hold off
    save(strcat(strrep(filename, '.csv', ''), '_class_Bakofen.mat'))
end