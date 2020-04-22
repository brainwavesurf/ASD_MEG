function [freq_with, freq_without, slope, intercept] = alpha_freq(x,y)
index = ((x>=2 & x<=5) | (x>=30 & x<=40));
coefficients = polyfit(log(x(index)), log(y(1,(index))), 1);
slope = coefficients(1,1);
intercept = coefficients(1,2);
trend = polyval(coefficients, log(x));

chis_with = 0; znam_with = 0;
chis_without = 0; znam_without = 0;

without = exp(log(y)) - exp(trend);

for f = 8:13
    chis_with = chis_with + (y(1,f))*f;
    znam_with = znam_with + (y(1,f));
    chis_without = chis_without + (without(1,f))*f;
    znam_without = znam_without + (without(1,f));
end 

freq_with = chis_with/znam_with; 
freq_without = chis_without/znam_without;
end
