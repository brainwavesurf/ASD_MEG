function freq = alpha_freq(x,y)
index = ((x>=2 & x<=5) | (x>=30 & x<=40));
coefficients = polyfit(log(x(index)), log(y(1,(index))), 1);
trend = polyval(coefficients, log(x));chis = 0;
znam = 0;
for f = 8:13
    without = exp(log(y)) - exp(trend);
    chis = chis + (without(1,f))*f;
    znam = znam + (without(1,f));
end  
freq = chis/znam;
end
