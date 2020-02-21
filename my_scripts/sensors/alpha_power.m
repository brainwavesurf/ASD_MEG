function power = alpha_power(x,y)
index = (x<=5 | x>=30);
coefficients = polyfit(log(x(index)), log(y(1,(index))), 1);
trend = polyval(coefficients, log(x));
for f = 8:13
    power(f) = max((log(y(1,f)) - trend));
end
power = sum(power)/6;
end