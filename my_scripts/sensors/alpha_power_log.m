function power = alpha_power_log(x,y)
index = ((x>=2 & x<=5) | (x>=30 & x<=40));
[slope, intercept] = logfit(x(index), y(1,index), 'logx');
trend = intercept + slope*log10(x);
for f = 8:13
    power(f) = y(1,f) - trend(1,f);
end
power = sum(power)/length(8:13);
end

        