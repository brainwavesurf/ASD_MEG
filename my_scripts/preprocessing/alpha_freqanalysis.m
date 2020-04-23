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
    load(strcat(savepath, subj, '/', subj, '_preproc_alpha_2_40_epochs.mat'));
    
    cfg = [];
    cfg.method         = 'mtmfft';
    cfg.output         = 'pow'; 
    cfg.taper          = 'hanning'; %Hanning taper
    %cfg.keeptrials     = 'yes';   
    cfg.pad            = 'nextpow2';
    cfg.foilim         = [2 40];
    
    fft_fast_isi        = ft_freqanalysis(cfg, fast_alpha_isi);
    fft_slow_isi        = ft_freqanalysis(cfg, slow_alpha_isi);
    fft_fast_post       = ft_freqanalysis(cfg, fast_alpha_post);
    fft_slow_post       = ft_freqanalysis(cfg, slow_alpha_post);
    
    filename = strcat(savepath, subj, '/', subj, '_alpha_freqanalysis.mat');
    save(filename, 'fft_fast_isi', 'fft_slow_isi', 'fft_fast_post', 'fft_slow_post');
end