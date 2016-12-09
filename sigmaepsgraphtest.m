A_mas_l = linspace(-100, 100, 10);
eps_dot = linspace(0, 0.001, 1000);
X_min = [19.524889977382376 1.812790936438664e+02 5.275591250583913 23.902573812023476];
for i = 1:10
    figure
    
    plot(0, 0)
    hold on
    for j = 1:10
        plot(eps_dot, A_mas_l(i)*log(eps_dot)+X_min(2) + def*(A_mas_l(j)*log(eps_dot)+X_min(4)), 'DisplayName', num2str(A_mas_l(j)));
        title(['A1 = ', num2str(A_mas_l(i))])
        xlabel('eps_dot')
        ylabel('sigma')
    end
    legend('show')
end