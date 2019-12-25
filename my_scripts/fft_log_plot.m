%plot log power after fft for all subjects individually
%in fast and slow conditions for grad and mag
%calculate coefficients for logarithmic fit and substract it from data

clear all;
close all;

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(path, fieldtripfolder)
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);
path('/home/a_shishkina/externals/', path);

savepath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';
realdatapath = '/home/a_shishkina/data/KI/SUBJECTS/';

%load subj list
     
SUBJ_NT = [ '0101'; '0102'; '0103'; '0104'; '0105'; '0135'; '0136';...  
            '0137'; '0138'; '0140'; '0158'; '0162'; '0163'; '0178';...
            '0179'; '0255'; '0257'; '0348'; '0378'; '0379'; '0384']; 
        
SUBJ_ASD = ['0106'; '0107'; '0139'; '0141'; '0159'; '0160'; '0161';...  
            '0164'; '0253'; '0254'; '0256'; '0273'; '0274'; '0275';...
            '0276'; '0346'; '0347'; '0349'; '0351'; '0358';...
            '0380'; '0381'; '0382'; '0383'];  
%without '0357';
SUBJ = [SUBJ_ASD; SUBJ_NT];

%load posterior sensors list
post_sens = {'MEG1932',  'MEG1922', 'MEG2042',  'MEG2032',  'MEG2112', 'MEG2122',  'MEG2342', 'MEG2332',  'MEG1732', 'MEG1942', 'MEG1912', 'MEG2012', 'MEG2022', 'MEG2312', 'MEG2322', 'MEG2512',...
             'MEG1933',  'MEG1923', 'MEG2043',  'MEG2033',  'MEG2113', 'MEG2123',  'MEG2343', 'MEG2333',  'MEG1733', 'MEG1943', 'MEG1913', 'MEG2013', 'MEG2023', 'MEG2313', 'MEG2323', 'MEG2513'};
    
%loop for all subjects
for s=1: size (SUBJ,1)
    
    subj = SUBJ (s,:); 
    savemegto = strcat(savepath, subj);
    epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
    
    %load preprocessed epochs
    epo = load(strcat(epofolder, subj, '_preproc_epochs.mat'));
    
    %select grad and mag epochs separately for slow and fast conditions
    cfg = [];
    cfg.channel = 'MEGMAG';
    epo_fast_mag = ft_selectdata(cfg, epo.fast_epochs);
    epo_slow_mag = ft_selectdata(cfg, epo.slow_epochs);
    
    cfg.channel = 'MEGGRAD';
    epo_fast_grad = ft_selectdata(cfg, epo.fast_epochs);
    epo_slow_grad = ft_selectdata(cfg, epo.slow_epochs);
    
 
    cfg = [];
    cfg.method       = 'mtmfft';
    cfg.output       = 'pow'; 
    cfg.taper        = 'hanning'; %Hanning taper
    cfg.foilim       = [5 30];       
    cfg.tapsmofrq    = 3; 
    fft_fast_grad{s} = ft_freqanalysis(cfg, epo_fast_grad);
    fft_slow_grad{s} = ft_freqanalysis(cfg, epo_slow_grad);
    fft_fast_mag{s} = ft_freqanalysis(cfg, epo_fast_mag);
    fft_slow_mag{s} = ft_freqanalysis(cfg, epo_slow_mag);
    
    %save stats
    filename = strcat(savepath, subj, '/', subj, '_fft_freq_analysis.mat');
    save(filename, 'fft_fast_grad', 'fft_slow_grad', 'fft_fast_mag', 'fft_slow_mag');
    
    %calculate coefficients for logarithmic fit 
    [slope1, intersep1] = logfit(fft_slow_grad{s}.freq, log(fft_slow_grad{s}.powspctrm(1,:)), 'logx');
    [slope2, intersep2] = logfit(fft_fast_grad{s}.freq, log(fft_fast_grad{s}.powspctrm(1,:)), 'logx');
    [slope3, intersep3] = logfit(fft_slow_mag{s}.freq, log(fft_slow_mag{s}.powspctrm(1,:)), 'logx');
    [slope4, intersep4] = logfit(fft_fast_mag{s}.freq, log(fft_fast_mag{s}.powspctrm(1,:)), 'logx');
    
    %plot figures for log power with linear trend 
    %for grad
    figure(1)
    subplot(2,1,1)
    plot(fft_slow_grad{s}.freq, log(fft_slow_grad{s}.powspctrm(1,:)), '-b'); title([subj, ', FFT power, grad'])
    hold on
    plot(intersep1 + slope1*log10(fft_fast_grad{s}.freq),'k--')
    legend('slow', 'lin trend'); xlim([5,30]); xlabel('Frequency (Hz)'); ylabel('log power');
    hold off 
    
    subplot(2,1,2)
    plot(fft_fast_grad{s}.freq, log(fft_fast_grad{s}.powspctrm(1,:)), '-r');
    hold on
    plot(intersep2 + slope2*log10(fft_slow_grad{s}.freq),'k--')
    legend('fast', 'lin trend'); xlim([5,30]); xlabel('Frequency (Hz)'); ylabel('log power');
    hold off 
    
    %for mag
    figure(2)
    subplot(2,1,1)
    plot(fft_slow_mag{s}.freq, log(fft_slow_mag{s}.powspctrm(1,:)), '-b'); title([subj, ', FFT power, mag'])
    hold on
    plot(intersep3 + slope3*log10(fft_fast_mag{s}.freq),'k--')
    legend('slow', 'lin trend'); xlim([5,30]); xlabel('Frequency (Hz)'); ylabel('log power');
    hold off 
    
    subplot(2,1,2)
    plot(fft_fast_mag{s}.freq, log(fft_fast_mag{s}.powspctrm(1,:)), '-r');
    hold on
    plot(intersep4 + slope4*log10(fft_slow_mag{s}.freq),'k--')
    legend('fast', 'lin trend'); xlim([5,30]); xlabel('Frequency (Hz)'); ylabel('log power');
    hold off 
    
    %name linear trend to use for plot 
    lintrend_slow_grad = intersep1 + slope1*log10(fft_slow_grad{s}.freq);
    lintrend_fast_grad = intersep2 + slope2*log10(fft_fast_grad{s}.freq);
    lintrend_slow_mag = intersep3 + slope3*log10(fft_slow_mag{s}.freq);
    lintrend_fast_mag = intersep4 + slope4*log10(fft_fast_mag{s}.freq);
    
    %plot power for fast and slow without linear trend
    %for grad
    figure(3)
    subplot(2,1,1)
    plot(fft_slow_grad{s}.freq, log(fft_slow_grad{s}.powspctrm(1,:)) - lintrend_slow_grad, '-b'); title([subj, ', log power with removal linear trend, grad'])
    hold on
    plot(fft_fast_grad{s}.freq, log(fft_fast_grad{s}.powspctrm(1,:)) - lintrend_fast_grad, '-r'); 
    legend('slow', 'fast'); xlim([5,30]); xlabel('Frequency (Hz)'); ylabel('log relative power');
    hold off
    %for mag
    subplot(2,1,2)
    plot(fft_slow_mag{s}.freq, log(fft_slow_mag{s}.powspctrm(1,:)) - lintrend_slow_mag, '-b'); title([subj, ', log power with removal linear trend, mag'])
    hold on
    plot(fft_fast_mag{s}.freq, log(fft_fast_mag{s}.powspctrm(1,:)) - lintrend_fast_mag, '-r'); 
    legend('slow', 'fast'); xlim([5,30]); xlabel('Frequency (Hz)'); ylabel('log relative power');
    hold off
    
    %save figures
    saveas(figure(1), [savepath, subj, '/', subj, '_FFT_log_power_grad.jpeg']);
    saveas(figure(2), [savepath, subj, '/', subj, '_FFT_log_power_mag.jpeg']);
    saveas(figure(3), [savepath, subj, '/', subj, '_FFT_log_relative_power.jpeg']);
 end
   
%save stats
filename = strcat(savepath, '1_results/', 'fft_freq_analysis.mat');
save(filename, 'fft_fast_grad', 'fft_slow_grad', 'fft_fast_mag', 'fft_slow_mag');


