% Load data epochs 
% Load info about the preceding trials

% Calculate time-freq representations by 
% 1) mtmconvol with Hanning taper with varying length of window
% 2) wavelet method with enhancement of the sensitivity to the non-phase locked activity
% 3) mtmfft with multitapers 
% 4) mtmfft with Hanning tapper

% Plot the time-frequency representation (TRF) for grad and mag in the each condition for
% 1) multiple channels (only grad)
% 2) for one channel with max power according to the previous plot (2043)
% 3) time-frequency topoplot for alpha range in interstimulus [-0.8 0] 

% Calculate and plot Fast-Slow power differences for grad and mag

%%
clear all;
close all;

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(fieldtripfolder, path);
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);
realdatapath = '/home/a_shishkina/data/KI/SUBJECTS/';
savepath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';

%%
%add list of subjects:
SUBJ = ['0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0136'; '0137'; '0138'; '0139'; '0140'; '0141'; '0158'; '0159'; '0160'; '0160'; '0161'; '0162'; '0163'; '0164'; '0178'; '0179'; '0253'; '0254'; '0255'; '0256'; '0257'; '0259'; '0273'; '0274'; '0275'; '0276'; '0277'; '0346'; '0347'; '0348'; '0350'; '0351'; '0357'; '0358'; '0378'; '0380'; '0381'; '0382'; '0383'; '0384']; 

%%
for s=1: size (SUBJ,1)
    close all
    subj = SUBJ (s,:); 
    savemegto = strcat(savepath, subj);
    epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
    
    %% Load unfiltered epochs. If you need to filter the data (e.g. for LCMV), import raw, not epochs.
    ep_fiff_file = strcat(epofolder, subj, '-noerror-lagcorrected-epo.fif')
    hdr = ft_read_header(ep_fiff_file);
    
    cfg = [];  
    cfg.dataset = ep_fiff_file;
    cfg.channel={'meg'};
    epochs = ft_preprocessing(cfg);
    
    %% load group preceding events in order to select the epochs later on
    load ([ savemegto, '/', subj, '_info.mat'])
    ev1 = find(allinfo.prev_stim_type==2); %epochs following slow
    ev2 = find(allinfo.prev_stim_type==4); %medium
    ev3 = find(allinfo.prev_stim_type==8); %fast
    prev = {ev1, ev2, ev3};
    
    for con = 1:3 % for conditions
        %select epochs according to preceding trials
        cfg = [];
        cfg.trials = prev{con}; 
        dataCon = ft_selectdata(cfg, epochs);

        %select interstimulus interval
        cfg = [];
        cfg.latency = [-0.8 0.0];
        dataPre = ft_selectdata(cfg, dataCon);

        %% do spectral analysis using Hanning taper, for freq < 30 Hz
        %==================================================================
        %Get PSD with Hanning taper

        cfg = [];
        cfg.method       = 'mtmfft';
        cfg.output       = 'pow'; 
        cfg.taper        = 'hanning'; %Hanning taper
        cfg.pad          = 2; % to [vertually] increase spectral resolution
        cfg.tapsmofrq    = 2;
        cfg.foi          = 10.5; %alpharange; 
        freqPre{con}     = ft_freqanalysis(cfg, dataPre);  %trials x Ch x freq
        powPre{con}      = freqPre{con}.powspctrm

    end
    
    diff_v3v1 = freqPre{1}; % use freqPre{1} structure, but change power for the 'difference'
    diff_v3v1.powspctrm = diff_v3v1.powspctrm';
    diff_v3v1.powspctrm = (powPre{3} - powPre{1})./(powPre{1} + powPre{3});
    powDiff = diff_v3v1.powspctrm;
    
    %% save stat info
    filename = strcat(savemegto, '/', subj, '_alpha_power.mat');
    save (filename, 'freqPre', 'powPre', 'powDiff');
    
end