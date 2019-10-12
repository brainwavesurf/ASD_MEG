%% Cluster-based permutation test for within-subjects experiments 
%(http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/#within-subjects-experiments)
%for magnetometers
%for gradiometers
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
SUBJ = ['0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0136'; '0137'; '0138'; '0139'; '0140'; '0141'; '0158'; '0159'; '0160'; '0161'; '0162'; '0163'; '0164'; '0178'; '0179'; '0253'; '0254'; '0255'; '0256'; '0257'; '0259'; '0273'; '0274'; '0275'; '0276'; '0277'; '0346'; '0347'; '0348'; '0350'; '0351'; '0357'; '0358'; '0378'; '0380'; '0381'; '0382'; '0383'; '0384']; 

for s=1: size (SUBJ,1)
    
    subj = SUBJ (s,:); 
    savemegto = strcat(savepath, subj);
    epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
    
    %load preprocessed epochs
    epo = load(strcat(epofolder, subj, '_preproc_epochs.mat'));
    
    %select mag epochs separately for slow and fast conditions
    cfg = [];
    cfg.channel = 'MEGMAG';
    epo_fast_mag = ft_selectdata(cfg, epo.fast_epochs);
    epo_slow_mag = ft_selectdata(cfg, epo.slow_epochs);
    
    %load data after freqanalysis 
    freq_all = load(strcat(savepath, subj, '/', subj, '_freqanalysis.mat'));
    %after mtmconvol
    freqFast_mag{s}  = freq_all.wvlts_fast_mag;
    freqSlow_mag{s}  = freq_all.wvlts_slow_mag;
end
%% do cluster-based permutation test

cfg = [];
cfg.channel          = {'MEGMAG'};
cfg.latency          = [-0.7 0];
cfg.frequency        = 10.5;
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
% specifies with which sensors other sensors can form clusters
cfg_neighb = [];
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, freqFast_mag{1}); %using a distance threshold of 4, there are on average 13.8 neighbours per channel

subj = size (SUBJ,1);
design = zeros(2,2*subj);
for i = 1:subj
design(1,i) = i;
end
for i = 1:subj
design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[within_subj_stat] = ft_freqstatistics(cfg, freqSlow_mag{:}, freqFast_mag{:})
stat_value = within_subj_stat.stat;

%save stats
filename = strcat(savepath, '1_results/', 'within_subjs_stat.mat');
save(filename, 'within_subj_stat');

filename = strcat(savepath, '1_results/', 'stat_value.mat');
save(filename, 'stat_value');
   
%plot
figure(1)
cfg = [];
cfg.alpha  = 0.025;
cfg.parameter = 'stat';
cfg.layout = 'neuromag306mag_helmet.mat';
ft_clusterplot(cfg, within_subj_stat);

%plot 
saveas(figure(1),[savepath, '/1_results/', 'cluster-based-permut-test.jpeg']);

%% the same for gradiometers

for s=1: size (SUBJ,1)
    
    subj = SUBJ (s,:); 
    savemegto = strcat(savepath, subj);
    epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
    
    %load preprocessed epochs
    epo = load(strcat(epofolder, subj, '_preproc_epochs.mat'));
    
    %select mag epochs separately for slow and fast conditions
    cfg = [];
    cfg.channel = 'MEGGRAD';
    epo_fast_grad = ft_selectdata(cfg, epo.fast_epochs);
    epo_slow_grad = ft_selectdata(cfg, epo.slow_epochs);
    
    %load data after freqanalysis 
    freq_all = load(strcat(savepath, subj, '/', subj, '_freqanalysis.mat'));
    %after mtmconvol
    freqFast_grad{s}  = freq_all.wvlts_fast_grad;
    freqSlow_grad{s}  = freq_all.wvlts_slow_grad;
end

%% do cluster-based permutation test
cfg = [];
cfg.channel          = {'MEGGRAD'};
cfg.latency          = [-0.7 0];
cfg.frequency        = 10.5;
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
% specifies with which sensors other sensors can form clusters
cfg_neighb = [];
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, freqFast_grad{1}); %using a distance threshold of 4, there are on average 13.8 neighbours per channel

subj = size (SUBJ,1);
design = zeros(2,2*subj);
for i = 1:subj
design(1,i) = i;
end
for i = 1:subj
design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[within_subj_stat_grad] = ft_freqstatistics(cfg, freqSlow_grad{:}, freqFast_grad{:})
stat_value = within_subj_stat.stat;

%save stats
filename = strcat(savepath, '1_results/', 'within_subjs_stat_grad.mat');
save(filename, 'within_subj_stat_grad');

filename = strcat(savepath, '1_results/', 'stat_value.mat');
save(filename, 'stat_value');
   
%plot
figure(2)
cfg = [];
cfg.alpha  = 0.025;
cfg.parameter = 'stat';
cfg.layout = 'neuromag306planar_helmet.mat';%for grad
ft_clusterplot(cfg, within_subj_stat_grad);

%plot 
saveas(figure(2),[savepath, '/1_results/', 'cluster_based_permut_test_grad.jpeg']);


