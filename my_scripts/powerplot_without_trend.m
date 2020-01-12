%plot power without linear trend for log-log case
%and without log trend for original data case

clear all;
close all;

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(path, fieldtripfolder)
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);
path('/home/a_shishkina/projects/ASD_MEG/Git/ASD_MEG/my_scripts', path);
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

SUBJ = ['0384'];
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
    cfg.latency = [-0.8 0.0];
    epo_fast_mag = ft_selectdata(cfg, epo.fast_epochs);
    epo_slow_mag = ft_selectdata(cfg, epo.slow_epochs);
    
    cfg.channel = 'MEGGRAD';
    epo_fast_grad = ft_selectdata(cfg, epo.fast_epochs);
    epo_slow_grad = ft_selectdata(cfg, epo.slow_epochs);
 
    cfg = [];
    cfg.method       = 'mtmfft';
    cfg.output       = 'pow'; 
    cfg.taper        = 'hanning'; 
    cfg.pad           = 2048/1000; 
    cfg.foilim       = [2 40];       
    cfg.tapsmofrq    = 3; 
    slow_avg_grad{s} = ft_freqanalysis(cfg, epo_fast_grad);
    fast_avg_grad{s} = ft_freqanalysis(cfg, epo_slow_grad);
    slow_avg_mag{s} = ft_freqanalysis(cfg, epo_fast_mag);
    fast_avg_mag{s} = ft_freqanalysis(cfg, epo_slow_mag);
    
    %save stats
    filename = strcat(savepath, subj, '/', subj, '_fft_freq_analysis.mat');
    save(filename, 'fft_fast_grad', 'fft_slow_grad', 'fft_fast_mag', 'fft_slow_mag');

    for i=1:2
        if i == 1
            type = 'loglog'; %log(freq) and log(power)
        else 
            type = 'original'; %freq and power
        end
        
        %for grad
        figure(1+3*(i-1))
        clf
        subplot_tight(2,1,1,[0.12 0.3])
        [trend_slow_grad] = plot_with_trend(slow_avg_grad{s}.freq,slow_avg_grad{s}.powspctrm,type,'grad','slow',subj);
        subplot_tight(2,1,2,[0.12 0.3])
        [trend_fast_grad] = plot_with_trend(fast_avg_grad{s}.freq,fast_avg_grad{s}.powspctrm,type,'grad','fast',subj);

        %for mag
        figure(1+3*(i-1)+1) 
        clf
        subplot_tight(2,1,1,[0.12 0.3])
        [trend_slow_mag] = plot_with_trend(slow_avg_mag{s}.freq,slow_avg_mag{s}.powspctrm,type,'mag','slow',subj);
        subplot_tight(2,1,2,[0.12 0.3])
        [trend_fast_mag] = plot_with_trend(fast_avg_mag{s}.freq,fast_avg_mag{s}.powspctrm,type,'mag','fast',subj);
        
        %plot power for fast and slow without trend
        
        figure(1+3*(i-1)+2)
        clf
        
        if strcmp(type,'loglog')
            %for grad
            subplot_tight(2,1,1,[0.12 0.3])
            plot(slow_avg_grad{s}.freq, log(slow_avg_grad{s}.powspctrm(1,:)) - trend_slow_grad, '-b'); title([subj, ', power with removal linear trend, grad'])
            hold on
            plot(fast_avg_grad{s}.freq, log(fast_avg_grad{s}.powspctrm(1,:)) - trend_fast_grad, '-r'); 
            legend('slow', 'fast'); xlabel('Frequency (Hz)'); ylabel('relative power');
            xlim([2,40])
            hold off
            set(gca,'Xscale','log')
            set(gca,'XTick',[5,10,20,30,40],...
                    'XTickLabel',{'5' '10' '20' '30' '40'})
            %for mag
            subplot_tight(2,1,2,[0.12 0.3])
            plot(slow_avg_mag{s}.freq, log(slow_avg_mag{s}.powspctrm(1,:)) - trend_slow_mag, '-b'); title([subj, ', power with removal linear trend, mag'])
            hold on
            plot(fast_avg_mag{s}.freq, log(fast_avg_mag{s}.powspctrm(1,:)) - trend_fast_mag, '-r'); 
            legend('slow', 'fast'); xlabel('Frequency (Hz)'); ylabel('relative power');
            xlim([2,40])
            hold off
            set(gca,'Xscale','log')
            set(gca,'XTick',[5,10,20,30,40],...
                'XTickLabel',{'5' '10' '20' '30' '40'})
        else
            %for grad
            subplot_tight(2,1,1,[0.12 0.3])
            plot(slow_avg_grad{s}.freq, slow_avg_grad{s}.powspctrm(1,:) - trend_slow_grad, '-b'); title([subj, ', power with removal log trend, grad'])
            hold on
            plot(fast_avg_grad{s}.freq, fast_avg_grad{s}.powspctrm(1,:) - trend_fast_grad, '-r'); 
            legend('slow', 'fast'); xlabel('Frequency (Hz)'); ylabel('relative power');
            xlim([2,40])
            hold off
             
            %for mag
            subplot_tight(2,1,2,[0.12 0.3])
            plot(slow_avg_mag{s}.freq, slow_avg_mag{s}.powspctrm(1,:) - trend_slow_mag, '-b'); title([subj, ', power with removal log trend, mag'])
            hold on
            plot(fast_avg_mag{s}.freq, fast_avg_mag{s}.powspctrm(1,:) - trend_fast_mag, '-r'); 
            legend('slow', 'fast'); xlabel('Frequency (Hz)'); ylabel('relative power');
            xlim([2,40])
            hold off
        end
    end
    
    %save figures
    saveas(figure(1), [savepath, subj, '/', subj, '_FFT_log_power_grad.jpeg']);
    saveas(figure(2), [savepath, subj, '/', subj, '_FFT_log_power_mag.jpeg']);
    saveas(figure(3), [savepath, subj, '/', subj, '_FFT_relative_log_power.jpeg']);
    saveas(figure(4), [savepath, subj, '/', subj, '_FFT_power_grad.jpeg']);
    saveas(figure(5), [savepath, subj, '/', subj, '_FFT_power_mag.jpeg']);
    saveas(figure(6), [savepath, subj, '/', subj, '_FFT_relative_power.jpeg']);
 end
   
%save stats
filename = strcat(savepath, '1_results/', 'fft_freq_analysis.mat');
save(filename, 'fft_fast_grad', 'fft_slow_grad', 'fft_fast_mag', 'fft_slow_mag');

