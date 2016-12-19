function [fun] = minim2(H_exp, H, R0, Rho0, t_exp, P)
    f2 = 0;
    f1 = 0;
    %delta = numel(H) - num0;
    delta = 0;
%     test = [];
%     for i = 1:numel(H_exp)
%         if (i<numel(H_exp) && i+3+delta<=numel(H) && i+delta+3<numel(P)) && (t_exp(i)>=275) && ((H_exp(i) < R0+Rho0) && (P(i-3+delta) == P(i+3+delta)))
% %             if i == 309
% %                 P(i-3+delta);
% %                 P(i+3+delta);
% %             end
%             slag = ((H_exp(i+3)-H_exp(i-3))/(t_exp(i+3) - t_exp(i-3))-(H(i+3 + delta)-H(i-3 + delta))/(t_exp(i+3) - t_exp(i-3)))^2;
%             f1 = f1 + ((H_exp(i+3)-H_exp(i-3))/(t_exp(i+3) - t_exp(i-3))-(H(i+3 + delta)-H(i-3 + delta))/(t_exp(i+3) - t_exp(i-3)))^2;
%             test = [test; slag];
%         else
%             test = [test; 0];
%         end
%         
%     end
%     f1 = sqrt(f1)*100;
     i = 1;
    while i<min(numel(H_exp), numel(H)) %&& (H(i) <= R0 || H_exp(i)<=R0)
        if (H_exp(i) < R0+Rho0) && (H_exp(i) ~= 0)
            f2 = f2 + abs(H_exp(i)-H(i+delta));
%             if abs(H_exp(i)-H(i+delta))>1;
%                 i
%             end
        end
        i = i + 1;
    end
%     fun = f1+f2;
fun = f2;
end

