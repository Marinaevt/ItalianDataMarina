function [ F ] = wisdom( x )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
F = 0;
N = 5000;
for i = 1:N
    F = F + (x(2) - x(1)).^2;
end
F = sqrt(F)/sqrt(N);
end

