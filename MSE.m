function [Err] = MSE(x, y)
n = numel(x);
Err = 0;
for i = 1:n
    Err = Err + (x(i) - y(i))^2;
end
    Err = sqrt(Err/n);
end

