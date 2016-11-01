function [fun] = minim3(H_exp, H, R0)
    fun = 0;
    i = 1;
        while i<numel(H) %&& (H(i) <= R0 || H_exp(i)<=R0)
            fun = fun + (H_exp(i)-H(i))^2;
            i = i + 1;
        end
    fun = sqrt(fun / i);
end

