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
    sref = 0.0000220769518118314;
    s_0 = 4.28166914057785;
    s_1 = 12.98367869;
    m0 = 0.654006407143157;
    mju = -0.54197793;
    R0 = 22.5;
    Rho0 = 3;
    
    P = a(:,1);
    P = round(P/0.25)*0.25;
    H_exp = a(:,2);
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
    for i = 1:numel(P)
        if P(i) ~= Pstart
            [T, Y] = ode45(@(t, h)Testdhdt5(P(j), s0, Rho0, R0, h, sf, sref, s_0, s_1, m0, mju), t(j:i), Y(end));
            p45 = plot(T, Y, 'm', 'DisplayName', 'Starting values');
            H2 = [H2; Y];
            hold on
            j = i;
            Pstart = P(i);
        end
    end
     [X_min, Err] = fminsearchbnd(@(x)minsearcher6(H_exp, P, s0, Rho0, R0, t, sf, x(1), x(2), x(3), x(4), x(5)), [sref, s_0, s_1, m0, mju], [inf inf inf inf inf], [inf inf inf 1 inf], optimset('Display', 'iter', 'MaxFunEvals', 3000));   
     X_min
     Err
     Pstart = P(1);
     Y(end) = 1;
     j = 1;
    H1 = [];
    for i = 1:numel(P)
        if P(i) ~= Pstart
            [T, Y] = ode45(@(t, h)Testdhdt5(P(j), s0, Rho0, R0, h, sf, X_min(1), X_min(2), X_min(3), X_min(4), X_min(5)), t(j:i), Y(end));
            pkmn = plot(T, Y, 'r', 'DisplayName', 'Inverse analysis values');
            hold on
            H1 = [H1; Y];
            j = i;
            Pstart = P(i);
        end
    end
    mineu = minim2(H_exp, h, R0, Rho0, t, P)
    minode = minim2(H_exp, H1, R0, Rho0, t, P)

    hold on
    pexp = plot(t, H_exp, 'g', 'DisplayName', 'Experimental data');
    xlabel('Time')
    ylabel('Height')
    legend([pexp p45 pkmn])
    axis([0 inf 0 35])
    savefig(strcat(strrep(filename, '.csv', ''), '_strain_dep_Bak.fig'))
    hold off
    save(strcat(strrep(filename, '.csv', ''), '_strain_dep_Bak.mat'))
end
