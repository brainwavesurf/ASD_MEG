
% Calculate time-freq representations by 
% mtmconvol with Hanning taper with varying length of window
% and spectral analysis with Fourier using Hanning taper

close all;

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(fieldtripfolder, path);
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);
realdatapath = '/net/server/data/Archive/aut_gamma/orekhova/KI/SUBJECTS/';
savepath = '/net/server/data/Archive/aut_gamma/orekhova/KI/Scripts_bkp/Shishkina/KI/Results_Alpha_and_Gamma/';

% load subj info  
SUBJ_NT = [ '0101'; '0102'; '0103'; '0104'; '0105'; '0135'; '0136';...  
            '0137'; '0138'; '0140'; '0158'; '0162'; '0163'; '0178';...
            '0179'; '0255'; '0257'; '0348'; '0378'; '0379'; '0384']; 
        
SUBJ_ASD = ['0106'; '0107'; '0139'; '0141'; '0159'; '0160'; '0161';...  
            '0164'; '0253'; '0254'; '0256'; '0273'; '0274'; '0275';...
            '0276'; '0346'; '0347'; '0349'; '0351'; '0358';...
            '0380'; '0381'; '0382'; '0383'];  

SUBJ = [SUBJ_ASD; SUBJ_NT];

%loop for all subjects
for s=1: size (SUBJ,1)
    
    subj = SUBJ (s,:); 
    %load preprocessed epochs bandpassed 2-40 Hz
    epo = load(strcat(savepath, subj, '/', subj, '_preproc_alpha_2_40_epochs.mat'));
    
    %select grad and mag epochs separately for slow and fast conditions
    
    cfg = [];
    cfg.channel = 'MEGMAG'; %magnetometers
    %whole epoch
    epo_fast_mag = ft_selectdata(cfg, epo.fast_alpha_epochs);
    epo_slow_mag = ft_selectdata(cfg, epo.slow_alpha_epochs);
    %interstimulus
    epo_fast_mag_pre = ft_selectdata(cfg, epo.fast_alpha_isi);
    epo_slow_mag_pre = ft_selectdata(cfg, epo.slow_alpha_isi);
    %stimulation
    epo_fast_mag_post = ft_selectdata(cfg, epo.fast_alpha_post);
    epo_slow_mag_post = ft_selectdata(cfg, epo.slow_alpha_post);
    
    cfg = [];
    cfg.channel = 'MEGGRAD'; %gradiometers
    %whole epoch
    epo_fast_grad = ft_selectdata(cfg, epo.fast_alpha_epochs);
    epo_slow_grad = ft_selectdata(cfg, epo.slow_alpha_epochs);
    %interstimulus
    epo_fast_grad_pre = ft_selectdata(cfg, epo.fast_alpha_isi);
    epo_slow_grad_pre = ft_selectdata(cfg, epo.slow_alpha_isi);
    %stimulation
    epo_fast_grad_post = ft_selectdata(cfg, epo.fast_alpha_post);
    epo_slow_grad_post = ft_selectdata(cfg, epo.slow_alpha_post);

    %% time-frequency response with Hann taper with varying length
    %mtmconvol gives the time-frequency representation of the trial
    cfg            = [];
    cfg.output     = 'pow';
    cfg.method     = 'mtmconvol';
    cfg.taper      = 'hanning';  
    cfg.foi        = 2:40;
    cfg.t_ftimwin  = 3./cfg.foi; 
    cfg.toi        = -0.8:0.05:1;
    cfg.pad        = 'nextpow2';
    conv_fast_mag  = ft_freqanalysis(cfg, epo_fast_mag);
    conv_slow_mag  = ft_freqanalysis(cfg, epo_slow_mag);
    conv_fast_grad = ft_freqanalysis(cfg, epo_fast_grad);
    conv_slow_grad = ft_freqanalysis(cfg, epo_slow_grad);
    
    %time-frequency plot of all channels for the one subject
    cfg = [];
    cfg.baseline     = [0.4 1.0];
    cfg.baselinetype = 'absolute';
    cfg.zlim         = 'maxabs';
    cfg.showlabels   = 'yes';
    cfg.layout       = 'neuromag306mag_helmet.mat';
    figure
    colorbar;
    ft_multiplotTFR(cfg, conv_fast_mag);
    saveas(figure(1),[savepath, '/1_results/freqanalysis/', subj, '_TFR_trial_all_chan_mag.jpeg']);
    
    %time-frequency plot of the one channel for the one subject
    cfg              = [];
    cfg.baseline     = [0.4 1.2];
    cfg.maskstyle    = 'saturation';
    cfg.zlim         = [-6e-26 6e-26];
    cfg.channel      = 'MEG1841';
    cfg.interactive  = 'no';
    cfg.layout       = 'neuromag306mag_helmet.mat';
    
    figure(2)
    subplot(1,2,1)
    ft_singleplotTFR(cfg, conv_slow_mag);
    title('"slow"')
    xlabel('time, s')
    ylabel('frequency, Hz')
    ylim([4 20])
    axis image; axis square; 
    
    subplot(1,2,2)
    ft_singleplotTFR(cfg, conv_fast_mag);
    title('"fast"')
    xlabel('time, s')
    ylabel('frequency, Hz')
    ylim([4 20])
    axis image; axis square; 
    saveas(figure(2),[savepath, '/1_results/freqanalysis/', subj, '_TFR_trial_one_chan_mag.jpeg']);
    
    %% spectral analysis with Fourier using Hanning taper

    cfg = [];
    cfg.method         = 'mtmfft';
    cfg.output         = 'pow'; 
    cfg.taper          = 'hanning'; %Hanning taper
    %cfg.keeptrials     = 'yes';   
    cfg.pad            = 'nextpow2';
    cfg.foilim         = [2 40];
    fft_fast_mag       = ft_freqanalysis(cfg, epo_fast_mag); 
    fft_slow_mag       = ft_freqanalysis(cfg, epo_slow_mag);
    fft_fast_grad      = ft_freqanalysis(cfg, epo_fast_grad);
    fft_slow_grad      = ft_freqanalysis(cfg, epo_slow_grad);
    
    fft_fast_mag_pre   = ft_freqanalysis(cfg, epo_fast_mag_pre); 
    fft_slow_mag_pre   = ft_freqanalysis(cfg, epo_slow_mag_pre);
    fft_fast_mag_post  = ft_freqanalysis(cfg, epo_fast_mag_post); 
    fft_slow_mag_post  = ft_freqanalysis(cfg, epo_slow_mag_post);
    
    fft_fast_grad_pre  = ft_freqanalysis(cfg, epo_fast_grad_pre);
    fft_slow_grad_pre  = ft_freqanalysis(cfg, epo_slow_grad_pre);
    fft_fast_grad_post = ft_freqanalysis(cfg, epo_fast_grad_post);
    fft_slow_grad_post = ft_freqanalysis(cfg, epo_slow_grad_post);
    
    cfg = [];
    cfg.avgoverchan = 'yes';
    avg_fft_fast_mag = ft_selectdata(cfg, fft_fast_mag);
    avg_fft_slow_mag = ft_selectdata(cfg, fft_slow_mag);
    
    figure(3);
    hold on;
    plot(avg_fft_fast_mag.freq, avg_fft_fast_mag.powspctrm)
    plot(avg_fft_slow_mag.freq, avg_fft_slow_mag.powspctrm)
    legend('"fast"','"slow"')
    xlabel('Frequency, Hz');
    ylabel('absolute power, uV^2');
    saveas(figure(3),[savepath, '/1_results/freqanalysis/', subj, '_Fourier_avg_chan_mag.jpeg']);
    
    %save freq analysis results
    filename = strcat(savepath, subj, '/', subj, '_freqanalysis.mat');
    save (filename, 'conv_fast_mag', 'conv_slow_mag', 'conv_fast_grad', 'fft_slow_grad', ...
        'fft_fast_mag', 'fft_slow_mag', 'fft_fast_grad', 'wvlts_slow_grad', ...
        'fft_fast_mag_pre', 'fft_slow_mag_pre', 'fft_fast_mag_post', 'fft_slow_mag_post', ...
        'fft_fast_grad_pre', 'fft_slow_grad_pre', 'fft_fast_grad_post', 'fft_slow_grad_post');

end

   