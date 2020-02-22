function power = alpha_power(x,y)
index = ((x>=2 & x<=5) | (x>=30 & x<=40));
coefficients = polyfit(log(x(index)), log(y(1,(index))), 1);
trend = polyval(coefficients, log(x));
for f = 8:13
    without = exp(log(y)) - exp(trend);
    power(f) = without(1,f);
end
power = sum(power)/length(8:13);
end
