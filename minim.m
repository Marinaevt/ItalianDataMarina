function [fun] = minim(H_exp, H)
    fun = ((H_exp-H)/H_exp)^2;
end

