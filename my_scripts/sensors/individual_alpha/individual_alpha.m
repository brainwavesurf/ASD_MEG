
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
 
    %do freqanalysis
    cfg = [];
    cfg.method       = 'mtmfft';
    cfg.output       = 'pow'; 
    cfg.taper        = 'hanning'; 
    cfg.pad           = 2048/epo_fast_grad.fsample; 
    cfg.foilim       = [2 40];       
    cfg.tapsmofrq    = 3; 
    slow_avg_grad{s} = ft_freqanalysis(cfg, epo_fast_grad);
    fast_avg_grad{s} = ft_freqanalysis(cfg, epo_slow_grad);
    slow_avg_mag{s} = ft_freqanalysis(cfg, epo_fast_mag);
    fast_avg_mag{s} = ft_freqanalysis(cfg, epo_slow_mag);
    
    % calculate individual alpha power
    alpha_slow_grad{s} = alpha_power(slow_avg_grad{s}.freq,slow_avg_grad{s}.powspctrm);
    alpha_fast_grad{s} = alpha_power(fast_avg_grad{s}.freq,fast_avg_grad{s}.powspctrm);
    alpha_slow_mag{s} = alpha_power(slow_avg_mag{s}.freq,slow_avg_mag{s}.powspctrm);
    alpha_fast_mag{s} = alpha_power(fast_avg_mag{s}.freq,fast_avg_mag{s}.powspctrm);
    
    %save 
    filename = strcat(savepath, subj, '/', subj, '_alpha_power.mat');
    save(filename, 'alpha_slow_grad', 'alpha_fast_grad', 'alpha_slow_mag', 'alpha_fast_mag');
end