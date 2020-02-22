function freq = alpha_freq_log(x,y)
index = ((x>=2 & x<=5) | (x>=30 & x<=40));
[slope, intercept] = logfit(x(index), y(1,index), 'logx');
trend = intercept + slope*log10(x);
chis = 0; znam = 0;
for f = 8:13
    chis = chis + (y(1,f) - trend(1,f))*f;
    znam = znam + (y(1,f) - trend(1,f));
end
freq = chis/znam;
end
   