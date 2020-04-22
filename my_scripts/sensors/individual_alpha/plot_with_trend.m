function [trend] = plot_with_trend(x,y,type,sensor,cond,subj)
index = ((x>=2 & x<=5) | (x>=30 & x<=40));
switch type; 
    case {'original'}
        titlename = [', FFT power, ', sensor];
        [slope, intercept] = logfit(x(index), y(1,index), 'logx');
        trend = intercept + slope*log10(x);
    case {'loglog'}
        titlename = [', log FFT power, ', sensor];
        coefficients = polyfit(log(x(index)), log(y(1,(index))), 1);
        trend = polyval(coefficients, log(x));
        y = log(y(1,:));
end

if cond == 'slow'
    color = '-b';
elseif cond == 'fast'
    color = '-r';
end

plot(x, y(1,:), color); title([subj, titlename])
hold on
plot(x, trend,'k--')
xlim([2,40])

switch type; 
    case {'original'}
        legend(cond, 'log trend'); xlabel('Frequency (Hz)'); ylabel('power');
    case {'loglog'}
        legend(cond, 'lin trend'); xlabel('Frequency (Hz)'); ylabel('log power');
        set(gca,'Xscale','log')
        set(gca,'XTick',[5,10,20,30,40],...
            'XTickLabel',{'5' '10' '20' '30' '40'})
end

hold off 
