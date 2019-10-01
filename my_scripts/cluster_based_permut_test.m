%% Cluster-based permutation test fo between-trial experiments (two cond)
% Reading the data for slow and fast conditions
% Calculation of time-frequency analysis
% Permutation test (posclusters: [1×2 struct], negclusters: [1×9 struct])
% Plotting the results      
% no clusters present with a p-value lower than the specified alpha, nothing to plot
%%
clear all;
close all;

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(path, fieldtripfolder)
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);

realdatapath = '/home/a_shishkina/data/KI/SUBJECTS/';
savepath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';

%%
%add list of subjects:
SUBJ = ['0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0136'; '0137'; '0138'; '0139'; '0140'; '0141'; '0158'; '0159'; '0160'; '0161'; '0162'; '0163'; '0164'; '0178'; '0179'; '0253'; '0254'; '0255'; '0256'; '0257'; '0259'; '0273'; '0274'; '0275'; '0276'; '0277'; '0346'; '0347'; '0348'; '0350'; '0351'; '0357'; '0358'; '0378'; '0380'; '0381'; '0382'; '0383'; '0384'; '0385']; 
%SUBJ = ['0076'];

for s=1: size (SUBJ,1)
    close all
    subj = SUBJ (s,:); 
    savemegto = strcat(savepath, subj);
    epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
    
    %% load group preceding events in order to select the epochs later on
    load ([ savemegto, '/', subj, '_info.mat'])
    
    slow_ind = find(allinfo.prev_stim_type==2);
    fast_ind = find(allinfo.prev_stim_type==8);
 
    %% Load unfiltered epochs. If you need to filter the data (e.g. for LCMV), import raw, not epochs.
    ep_fiff_file = strcat(epofolder, subj, '-noerror-lagcorrected-epo.fif')
    hdr = ft_read_header(ep_fiff_file);
    
    cfg           = [];  
    cfg.dataset   = ep_fiff_file;
    cfg.channel   = {'MEGGRAD'};
    epochs        = ft_preprocessing(cfg);
    
    cfg           = [];
    cfg.trials    = slow_ind;
    fast_epochs   = ft_selectdata(cfg, epochs);
    
    cfg.trials    = fast_ind;
    slow_epochs   = ft_selectdata(cfg, epochs);
    
    %% Calculate the TFRs for the two experimental conditions (Fast and Slow).
    cfg = [];
    cfg.output     = 'pow';
    cfg.channel    = 'MEG';
    cfg.method     = 'mtmconvol';
    cfg.taper      = 'hanning';
    cfg.foi        = [7:20];
    cfg.toi        = [-0.8:0.05:0];
    cfg.t_ftimwin  = 7./cfg.foi; %7 cycles
    cfg.keeptrials = 'yes';

    freqFast  = ft_freqanalysis(cfg, fast_epochs);
    freqSlow  = ft_freqanalysis(cfg, slow_epochs);

    %% do permutation test
    cfg = [];
    cfg.channel          = {'MEG'};
    cfg.latency          = 'all';
    cfg.frequency        = [7 14];
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'ft_statfun_indepsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan        = 2;
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.alpha            = 0.025;
    cfg.numrandomization = 500;
    % prepare_neighbours determines what sensors may form clusters
    cfg_neighb.method    = 'distance';
    cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, fast_epochs);

    design = zeros(1,size(freqSlow.powspctrm,1) + size(freqFast.powspctrm,1));
    design(1,1:size(freqSlow.powspctrm,1)) = 1;
    design(1,(size(freqSlow.powspctrm,1)+1):(size(freqSlow.powspctrm,1)+...
    size(freqFast.powspctrm,1))) = 2;

    cfg.design           = design;
    cfg.ivar             = 1;

    [stat] = ft_freqstatistics(cfg, freqSlow, freqFast);
    
    filename = strcat(savemegto, '/', subj, '_cluster_based_permut_sensors.mat');
    save (filename, 'stat');
    
    % Averaging over epochs
    cfg = [];
    cfg.avgoverfreq = 'yes';
    stat_avg = ft_selectdata(cfg, stat);
    
    %% Plotting the results
    
    cfg = [];
    freqSlow = ft_freqdescriptives(cfg, freqSlow);
    freqFast = ft_freqdescriptives(cfg, freqFast);
    
    cfg = [];
    cfg.avgoverfreq = 'yes';
    freqSlow_avg = ft_selectdata(cfg, freqSlow);
    freqFast_avg = ft_selectdata(cfg, freqFast);
    
    stat_avg.raweffect = freqFast_avg.powspctrm - freqSlow_avg.powspctrm;
    
    cfg = [];
    cfg.alpha  = 0.025;
    cfg.parameter = 'raweffect';
    cfg.zlim   = [-1e-27 1e-27];
    cfg.layout = 'neuromag306planar_helmet.mat';
    ft_clusterplot(cfg, stat_avg);

end
