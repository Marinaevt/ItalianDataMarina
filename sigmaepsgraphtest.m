A_mas_l = linspace(-100, 100, 10);
for i = 1:10
    figure
    title('A1 = %f', A_mas_l)
    xlabel('eps_dot')
    ylabel('sigma')
    plot(0, 0)
    hold on
    for j = 1:10
        plot(eps_dot, A_mas_l(i)*log(eps_dot)+X_min(2) + def*(A_mas_l(j)*log(eps_dot)+X_min(4)));
    end
end