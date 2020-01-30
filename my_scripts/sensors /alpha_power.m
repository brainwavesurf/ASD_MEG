function [alpha] = alpha_power(x,y)
index = (x<=5 | x>=30);
coefficients = polyfit(log(x(index)), log(y(1,(index))), 1);
trend = polyval(coefficients, log(x));chis = 0;
znam = 0;
for f = 8:13
    chis = chis + ((log(y(1,f)) - trend)*f);
    znam = znam + (log(y(1,f)) - trend);
end  
alpha = chis/znam;
end
  
 