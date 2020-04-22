function [power_with, power_without] = alpha_power(x,y)
index = ((x>=2 & x<=5) | (x>=30 & x<=40));
coefficients = polyfit(log(x(index)), log(y(1,(index))), 1);
trend = polyval(coefficients, log(x));
without = exp(log(y)) - exp(trend);
for f = 10:17
    power_with = y(1,f);
    power_without(f) = without(1,f);
end
power_with = sum(power_with)*1e25/length(10:17);
power_without = sum(power_without)*1e25/length(10:17);
end
