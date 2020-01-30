
% Calculate time-freq representations by 
% 1) mtmconvol with Hanning taper with varying length of window
% 2) wavelet analysis
% 3) mtmfft with Hanning tapper

%%
clear all;
close all;

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(fieldtripfolder, path);
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);
realdatapath = '/home/a_shishkina/data/KI/SUBJECTS/';
savepath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';

%% load subj info
     
SUBJ_NT = [ '0101'; '0102'; '0103'; '0104'; '0105'; '0135'; '0136';...  
            '0137'; '0138'; '0140'; '0158'; '0162'; '0163'; '0178';...
            '0179'; '0255'; '0257'; '0348'; '0378'; '0379'; '0384']; 
        
SUBJ_ASD = ['0106'; '0107'; '0139'; '0141'; '0159'; '0160'; '0161';...  
            '0164'; '0253'; '0254'; '0256'; '0273'; '0274'; '0275';...
            '0276'; '0346'; '0347'; '0349'; '0351'; '0358';...
            '0380'; '0381'; '0382'; '0383'];  
%without '0357';
SUBJ = [SUBJ_ASD; SUBJ_NT];
%% loop for all subjects
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
    
    cfg.latency = [-0.8 0];
    epo_fast_mag_pre = ft_selectdata(cfg, epo.fast_epochs);
    epo_slow_mag_pre = ft_selectdata(cfg, epo.slow_epochs);
    
    cfg.latency = [0.4 1.2];
    epo_fast_mag_post = ft_selectdata(cfg, epo.fast_epochs);
    epo_slow_mag_post = ft_selectdata(cfg, epo.slow_epochs);
    
    cfg = [];
    cfg.channel = 'MEGGRAD';
    epo_fast_grad = ft_selectdata(cfg, epo.fast_epochs);
    epo_slow_grad = ft_selectdata(cfg, epo.slow_epochs);
    
    cfg.latency = [-0.8 0];
    epo_fast_grad_pre = ft_selectdata(cfg, epo.fast_epochs);
    epo_slow_grad_pre = ft_selectdata(cfg, epo.slow_epochs);
    
    cfg.latency = [0.4 1.2];
    epo_fast_grad_post = ft_selectdata(cfg, epo.fast_epochs);
    epo_slow_grad_post = ft_selectdata(cfg, epo.slow_epochs);
      
    %% time-frequency response with Hann taper with varying length
    %mtmconvol gives the time-frequency representation of the trial

    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.foi          = [5:30];                 
    cfg.toi          = [-1:0.05:1.2]; 
    cfg.t_ftimwin    = 3./cfg.foi; 
    conv_fast_mag  = ft_freqanalysis(cfg, epo_fast_mag);
    conv_slow_mag  = ft_freqanalysis(cfg, epo_slow_mag);
    conv_fast_grad = ft_freqanalysis(cfg, epo_fast_grad);
    conv_slow_grad = ft_freqanalysis(cfg, epo_slow_grad);
    
    %% TFR with Morlet wavelets 

    cfg = [];
    cfg.method       = 'wavelet';
    cfg.foilim       = [5 30];                 
    cfg.toi          = [-1:0.05:1.2];
    cfg.width        = 3;
  
    wvlts_fast_mag  = ft_freqanalysis(cfg, epo_fast_mag);
    wvlts_slow_mag  = ft_freqanalysis(cfg, epo_slow_mag);
    wvlts_fast_grad = ft_freqanalysis(cfg, epo_fast_grad);
    wvlts_slow_grad = ft_freqanalysis(cfg, epo_slow_grad);
        
        
    %% spectral analysis using Hanning taper, for freq < 30 Hz

    cfg = [];
    cfg.method         = 'mtmfft';
    cfg.output         = 'pow'; 
    cfg.taper          = 'hanning'; %Hanning taper
    cfg.keeptrials     = 'yes';          
    cfg.pad            = 10; 
    cfg.foilim         = [5 30];
    fft_fast_mag_pre   = ft_freqanalysis(cfg, epo_fast_mag_pre); % trials x Ch x freq
    fft_slow_mag_pre   = ft_freqanalysis(cfg, epo_slow_mag_pre);
    fft_fast_mag_post  = ft_freqanalysis(cfg, epo_fast_mag_post); % trials x Ch x freq
    fft_slow_mag_post  = ft_freqanalysis(cfg, epo_slow_mag_post);
    
    fft_fast_grad_pre  = ft_freqanalysis(cfg, epo_fast_grad_pre);
    fft_slow_grad_pre  = ft_freqanalysis(cfg, epo_slow_grad_pre);
    fft_fast_grad_post = ft_freqanalysis(cfg, epo_fast_grad_post);
    fft_slow_grad_post = ft_freqanalysis(cfg, epo_slow_grad_post);
    
    %save freq analysis results
    filename = strcat(savepath, subj, '/', subj, '_freqanalysis.mat');
    save (filename, 'conv_fast_mag', 'conv_slow_mag', 'conv_fast_grad', 'conv_slow_grad', ...
        'wvlts_fast_mag', 'wvlts_slow_mag', 'wvlts_fast_grad', 'wvlts_slow_grad', ...
        'fft_fast_mag_pre', 'fft_slow_mag_pre', 'fft_fast_mag_post', 'fft_slow_mag_post', ...
        'fft_fast_grad_pre', 'fft_slow_grad_pre', 'fft_fast_grad_post', 'fft_slow_grad_post');

end

   