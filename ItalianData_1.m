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
    H_exp = a(:,3);
    t_exp = a(:,4);
    F = [];
    F1 = [];
    num0 = find(H_exp == 0);
    num0 = num0(1) - 1;
    num10 = find(t == 10)-1;
    pit = P(num10+1);
    k = 1;
    wisdom1 = @(x)(((x(2) - x(1))/x(2)).^2 + ((x(4) - x(3))/x(4)).^2); 
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
    plot(t, h, 'b', 'DisplayName', 'Calculated data')
    hold on
%      plot(F(:, 3), F(:, 1), 'r', 'DisplayName', 'Data considering error')
    plot(t_exp(1:num0), H_exp(1:num0), 'g', 'DisplayName', 'Experimental data');
    err = MSE(H_exp(1:num0), h(num10:numel(h)));
    err(1:numel(t)) = err;
    plot(t, err, 'r', 'DisplayName', 'MSE');
    err(1:numel(t)) = err.^2;
    plot(t, err, 'y', 'DisplayName', 'DX');
%     plot(t_exp(1:num0), F1(:, 1), 'y', 'DisplayName', 'Data considering error ^2')
%    plot(t_exp, H_exp, 'g', 'DisplayName', 'Experimental data');
    xlabel('Time')
    ylabel('Height')
    legend('show')
    hold off
    axis([0 inf 0 35])
    figure('Name', char(name))
    plot(t, SR)
    xlabel('Time')
    ylabel('SR')
    axis([0 inf 0 0.01])
end
