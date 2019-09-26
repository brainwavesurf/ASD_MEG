%Preprocessing the response-locked data
%Single-trial time-frequency responses
%Log-transform the single-trial power
%Compute the neighbours
%Compute the statistics
%Visualise the results
%%
clear all;
close all;

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(fieldtripfolder, path);
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);
realdatapath = '/home/a_shishkina/data/KI/SUBJECTS/';
savepath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';

%SUBJ = ['0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0136'; '0137'; '0138'; '0139'; '0140'; '0141'; '0158'; '0159'; '0160'; '0161'; '0162'; '0163'; '0164'; '0178'; '0179'; '0253'; '0254'; '0255'; '0256'; '0257'; '0259'; '0273'; '0274'; '0275'; '0276'; '0277'; '0346';'0347';'0348'; '0350'; '0351'; '0357'; '0358'; '0378'; '0380'; '0381'; '0382'; '0383'; '0384']; 
SUBJ = ['0076'];

for s=1: size (SUBJ,1)
    close all
    subj = SUBJ (s,:); 
    savemegto = strcat(savepath, subj);
    epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
    
    %% Load unfiltered epochs. If you need to filter the data (e.g. for LCMV), import raw, not epochs.
    ep_fiff_file = strcat(epofolder, subj, '-noerror-lagcorrected-epo.fif')
    hdr = ft_read_header(ep_fiff_file);
    
    cfg                     = [];  
    cfg.dataset             = ep_fiff_file;
    cfg.trials              = [find(allinfo.prev_stim_type==2), find(allinfo.prev_stim_type==8)];
    cfg.channel             ={'meg'};
    epochs                  = ft_preprocessing(cfg);

    %% load group preceding events in order to select the epochs later on
    load ([ savemegto, '/', subj, '_info.mat'])
    
    %              event: [1×227 int32]
    %            epo_num: [1×227 int64]
    %                 rt: [1×227 double]
    %     prev_stim_type: [1×227 int64]
    %   prev_stim_length: [1×227 int64]
    %time_to_previous_ev: [1×227 double]

    %% For statistical testing we need to compute single-trial time-frequency responses 
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.trials       = [find(allinfo.prev_stim_type==2), find(allinfo.prev_stim_type==8)];
    cfg.toi          = -0.8 : 0.01 : 0.0;
    cfg.foi          = 7:20;
    cfg.t_ftimwin    = ones(size(cfg.foi)) * 0.5;
    cfg.keeptrials   = 'yes';  % keep the TFR on individual trials, dimord: 'rpt_chan_freq_time'
    TFR_all          = ft_freqanalysis(cfg, epochs);
    
    %% Log-transform the single-trial power
    %Spectral power is not normally distributed. 
    %Although this is in theory not a problem for the non-parametric statistical test, 
    %its sensitivity is usually increased by log-transforming the values in the power spectrum.    
    cfg           = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'log10';
    TFR_logpow    = ft_math(cfg, TFR_all);

    %% to form the clusters we have to construct an explicit description of the neighbourhood of each channel  
    cfg           = [];
    cfg.channel   = 'MEG*1';
    cfg.method    = 'triangulation';
    cfg.grad      = TFR_all.grad;
    cfg.feedback  = 'yes';
    neighbours    = ft_prepare_neighbours(cfg);
    
    %% Compute the statistics without MCP correction  
    cfg           = [];
    cfg.channel   = 'MEG*1';
    cfg.statistic = 'indepsamplesT';
    cfg.ivar      = 1;
    cfg.design    = zeros(1, size(TFR_all.cfg.trials,2));

    cfg.design(1:51) = 1; %slow
    cfg.design(52:end) = 2; %fast

    cfg.method    = 'analytic';
    cfg.correctm  = 'no';
    TFR_stat     = ft_freqstatistics(cfg, TFR_logpow);

    %% do correction for multiple comparisons by three methods 
    cfg.method    = 'analytic';
    cfg.correctm  = 'bonferoni';
    TFR_stat_bonferoni     = ft_freqstatistics(cfg, TFR_logpow);

    cfg.method    = 'analytic';
    cfg.correctm  = 'fdr';
    TFR_stat_fdr     = ft_freqstatistics(cfg, TFR_logpow);

    cfg.method            = 'montecarlo';
    cfg.correctm          = 'cluster';
    cfg.numrandomization  = 1000; 
    cfg.neighbours        = neighbours;
    TFR_stat_montecarlo     = ft_freqstatistics(cfg, TFR_logpow);
    
    %Plot the results
    cfg               = [];
    cfg.marker        = 'on';
    cfg.layout        = 'neuromag306mag.lay';
    cfg.channel       = 'MEG*1';
    cfg.parameter     = 'stat';  % plot the t-value
    cfg.maskparameter = 'mask';  % use the thresholded probability to mask the data
    cfg.maskstyle     = 'saturation';

    figure(1); ft_multiplotTFR(cfg, TFR_stat);
    title('Stat without MCP correction'); colorbar;
    saveas(figure(1),[savepath, subj, '/', subj, '_stat_without_MCP_correction.jpeg' ]);
    
    figure(2); ft_multiplotTFR(cfg, TFR_stat_bonferoni);
    title('Stat with bonferoni correction'); colorbar;
    saveas(figure(2),[savepath, subj, '/', subj, '_stat_with_bonnferoni_correction.jpeg' ]);
    
    figure(3); ft_multiplotTFR(cfg, TFR_stat_fdr);
    title('Stat with fdr correction'); colorbar;
    saveas(figure(3),[savepath, subj, '/', subj, '_stat_with_fdr_correction.jpeg' ]);
    
    figure(4); ft_multiplotTFR(cfg, TFR_stat_montecarlo);
    title('Stat with montecarlo correction'); colorbar;
    saveas(figure(4),[savepath, subj, '/', subj, '_stat_with_montecarlo_correction.jpeg' ]);
    
    %% save stat info
    filename = strcat(savemegto, '/', subj, '_stat_sensors.mat');
    save (filename, 'TFR_stat', 'TFR_stat_bonferoni', 'TFR_stat_fdr', 'TFR_stat_montecarlo');

    
end