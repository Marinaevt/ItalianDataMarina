function [fun] = minim1(H_exp, H)
    fun = abs((H_exp-H)/H_exp);
end

