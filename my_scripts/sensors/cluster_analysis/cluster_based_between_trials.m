%% Cluster-based permutation test for between-trial experiments 
%(http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/#within-subjects-experiments)
% my questions are commented and boundared by "===", also I added info about stuctures 
% I have got results after permutation test stats (9 posclusters, 5 negclusters), but I cannot vizualize it by ft_clusterplot:
%Error using ft_clusterplot (line 208)
%no clusters present with a p-value lower than the specified alpha, nothing to plot
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
%SUBJ = ['0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0136'; '0137'; '0138'; '0139'; '0140'; '0141'; '0158'; '0159'; '0160'; '0161'; '0162'; '0163'; '0164'; '0178'; '0179'; '0253'; '0254'; '0255'; '0256'; '0257'; '0259'; '0273'; '0274'; '0275'; '0276'; '0277'; '0346'; '0347'; '0348'; '0350'; '0351'; '0357'; '0358'; '0378'; '0380'; '0381'; '0382'; '0383'; '0384'; '0385']; 
SUBJ = ['0102'];

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
    
    freq_all = load(strcat(savepath, subj, '/', subj, '_freqanalysis.mat'));
    
    freqFast_mag  = freq_all.conv_fast_mag;
    freqSlow_mag  = freq_all.conv_slow_mag;
    
    %% Do between trials stats
    cfg = [];
    cfg.channel          = {'MEGMAG'};
    cfg.latency          = [-0.8 0];
    cfg.frequency        = [5 20];
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
    
    cfg_neighb.method    = 'distance';
    cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, epo_fast_mag);

    design = zeros(1,size(freqSlow_mag.powspctrm,1) + size(freqFast_mag.powspctrm,1));
    design(1,1:size(freqSlow_mag.powspctrm,1)) = 1;
    design(1,(size(freqSlow_mag.powspctrm,1)+1):(size(freqSlow_mag.powspctrm,1)+...
    size(freqFast_mag.powspctrm,1))) = 2;
    
    cfg.design           = design;
    cfg.ivar             = 1;
    
    [stat] = ft_freqstatistics(cfg, freqSlow_mag, freqFast_mag);

% stat = 
% 
%   struct with fields:
% 
%                    prob: [102×16×17 double]
%             posclusters: [1×9 struct]
%     posclusterslabelmat: [102×16×17 double]
%         posdistribution: [1×500 double]
%             negclusters: [1×5 struct]
%     negclusterslabelmat: [102×16×17 double]
%         negdistribution: [1×500 double]
%                 cirange: [102×16×17 double]
%                    mask: [102×16×17 logical]
%                    stat: [102×16×17 double]
%                     ref: [102×16×17 double]
%                  dimord: 'chan_freq_time'
%                    freq: [1×16 double]
%                    grad: [1×1 struct]
%                   label: {102×1 cell}
%                    time: [1×17 double]
%                     cfg: [1×1 struct]

    %% do average over trials
    
    cfg = [];
    cfg.latency = [-0.8 0];
    freqFast_prestim = ft_selectdata(cfg, freqFast_mag);
    freqSlow_prestim = ft_selectdata(cfg, freqSlow_mag);
      
% freqFast_prestim = struct with fields:
%         label: {102×1 cell}
%          freq: [1×16 double]
%          time: [1×17 double]
%     powspctrm: [63×102×16×17 double]
%     cumtapcnt: [63×16 double]
%          grad: [1×1 struct]
%           cfg: [1×1 struct]
%        dimord: 'rpt_chan_freq_time'
    
    cfg = [];
    cfg.avgoverrpt = 'yes';
    freqFast_avg = ft_selectdata(cfg, freqFast_prestim);
    freqSlow_avg = ft_selectdata(cfg, freqSlow_prestim);
    
% freqFast_avg = struct with fields:
%         label: {102×1 cell}
%          freq: [1×16 double]
%          time: [1×17 double]
%     powspctrm: [102×16×17 double]
%          grad: [1×1 struct]
%           cfg: [1×1 struct]
%        dimord: 'chan_freq_time'

    %% add the raw effect (Fast-Slow) to the obtained stat structure and plot the largest cluster overlayed on the raw effect.
    stat.raweffect = freqFast_avg.powspctrm - freqSlow_avg.powspctrm; %raweffect: [102×16×17 double]
    
    cfg = [];
    cfg.alpha  = 0.025;
    cfg.parameter = 'raweffect';
    cfg.layout = 'neuromag306mag_helmet.mat';
    ft_clusterplot(cfg, stat);
    
    %Error using ft_clusterplot (line 155)
    %this only works if either frequency or time is a singleton dimension
    %======================================================================
    %to solve this problem a do averaging for frequencies. Is it correctly?
    %I used advice from https://mailman.science.ru.nl/pipermail/fieldtrip/2016-May/010415.html
    %Are there other ways to fix it?
    %======================================================================
    
    cfg  = [];
    cfg.avgoverfreq = 'yes';
    stat_avg_freq = ft_selectdata(cfg, stat);

% stat_avg_freq = 
%   struct with fields:
% 
%                    prob: [102×1×17 double]
%             posclusters: [1×9 struct]
%     posclusterslabelmat: [102×1×17 double]
%         posdistribution: [1×500 double]
%             negclusters: [1×5 struct]
%     negclusterslabelmat: [102×1×17 double]
%         negdistribution: [1×500 double]
%                 cirange: [102×1×17 double]
%                    mask: [102×1×17 double]
%                    stat: [102×1×17 double]
%                     ref: [102×1×17 double]
%                    freq: 12.4896
%                    grad: [1×1 struct]
%                   label: {102×1 cell}
%                    time: [1×17 double]
%                     cfg: [1×1 struct]
%               raweffect: [102×1×17 double]
%                  dimord: 'chan_freq_time'

    cfg = [];
    cfg.alpha  = 0.025;
    cfg.parameter = 'raweffect';
    cfg.layout = 'neuromag306mag_helmet.mat';
    ft_clusterplot(cfg, stat_avg_freq);

    %Error using ft_clusterplot (line 208)
    %no clusters present with a p-value lower than the specified alpha, nothing to plot
    
    %save stats
    filename = strcat(savepath, '1_results/', '_cluster_based_between_trails.mat');
    save(filename, 'stat', 'stat_avg_freq');
end 
