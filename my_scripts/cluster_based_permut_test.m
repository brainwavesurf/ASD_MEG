%% Cluster-based permutation test for within-subjects experiments 
%(http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/#within-subjects-experiments)
%for magnetometers
%for gradiometers
%control, ASD and combined
%[-0.7 0] to avoid NAN stat
% avgoverfreq alpha [8 12]
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
% SUBJ = [ '0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0136'; '0137'; '0138'; '0139';...  
%          '0140'; '0141'; '0158'; '0159'; '0160'; '0161'; '0162'; '0163'; '0164'; '0178'; '0179'; '0253';...
%          '0254'; '0255'; '0256'; '0257'; '0259'; '0273'; '0274'; '0275'; '0276'; '0346'; '0347'; '0348';... 
%          '0351'; '0357'; '0358'; '0378'; '0380'; '0381'; '0382'; '0383'; '0384']; 
 
     % exclude 0310 (female, researcher, no MRI), 0076(researcher), 0350(researcher/ASD), 0135 (bad
     % coregistration), 0277 (Female, researcher), 0278/0279 (mom/son-with-epi 'controls'), 349 (ASD/no MRI),
     % 0379 (NT no MRI), 
SUBJ_NT = [ '0101'; '0102'; '0103'; '0104'; '0105'; '0136';...  
            '0137'; '0138'; '0140'; '0158'; '0162'; '0163'; '0178';...
            '0179'; '0255'; '0257'; '0348'; '0378'; '0384']; 
        
SUBJ_ASD = ['0106'; '0107'; '0139'; '0141'; '0159'; '0160'; '0161';...  
            '0164'; '0253'; '0254'; '0256'; '0273'; '0274'; '0275';...
            '0276'; '0346'; '0347'; '0351'; '0357'; '0358';...
            '0380'; '0381'; '0382'; '0383'];  
        
SUBJ_ALL = [SUBJ_ASD; SUBJ_NT];
%% load freqanalysis data for all groups
for s=1:size(SUBJ_ALL,1)
    
    subj = SUBJ_ALL(s,:); 
    %load data after freqanalysis 
    freq_all = load(strcat(savepath, subj, '/', subj, '_freqanalysis.mat'));
    %after mtmconvol
    freqFast_mag{s}  = freq_all.wvlts_fast_mag;
    freqSlow_mag{s}  = freq_all.wvlts_slow_mag;
end

%% load freqanalysis data for control group
for s=1:size(SUBJ_NT,1)
    
    subj = SUBJ_NT(s,:); 
    %load data after freqanalysis 
    freq_all = load(strcat(savepath, subj, '/', subj, '_freqanalysis.mat'));
    %after mtmconvol
    freqFast_mag_NT{s}  = freq_all.wvlts_fast_mag;
    freqSlow_mag_NT{s}  = freq_all.wvlts_slow_mag;
end

%% load freqanalysis data for asd group
for s=1:size(SUBJ_ASD,1)
    
    subj = SUBJ_ASD(s,:); 
    %load data after freqanalysis 
    freq_all = load(strcat(savepath, subj, '/', subj, '_freqanalysis.mat'));
    %after mtmconvol
    freqFast_mag_ASD{s}  = freq_all.wvlts_fast_mag;
    freqSlow_mag_ASD{s}  = freq_all.wvlts_slow_mag;
end
%% do cluster-based permutation test for ALL

cfg = [];
cfg.channel          = {'MEGMAG'};
cfg.latency          = [-0.7 0];
cfg.frequency        = [8 12];
cfg.avgoverfreq      = 'yes';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
% specifies with which sensors other sensors can form clusters
cfg_neighb = [];
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, freqFast_mag{1}); %using a distance threshold of 4, there are on average 13.8 neighbours per channel

subj = size(SUBJ_ALL,1);
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

[within_subj_stat] = ft_freqstatistics(cfg, freqFast_mag{:}, freqSlow_mag{:})

%% do cluster-based permutation test for NT

cfg = [];
cfg.channel          = {'MEGMAG'};
cfg.latency          = [-0.7 0];
cfg.frequency        = [8 12];
cfg.avgoverfreq      = 'yes';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
% specifies with which sensors other sensors can form clusters
cfg_neighb = [];
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, freqFast_mag_NT{1}); %using a distance threshold of 4, there are on average 13.8 neighbours per channel

subj = size(SUBJ_NT,1);
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

[within_subj_stat_NT] = ft_freqstatistics(cfg, freqFast_mag_NT{:}, freqSlow_mag_NT{:})

%% do cluster-based permutation test for ASD

cfg = [];
cfg.channel          = {'MEGMAG'};
cfg.latency          = [-0.7 0];
cfg.frequency        = [8 12];
cfg.avgoverfreq      = 'yes';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
% specifies with which sensors other sensors can form clusters
cfg_neighb = [];
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, freqFast_mag_ASD{1}); %using a distance threshold of 4, there are on average 13.8 neighbours per channel

subj = size(SUBJ_ASD,1);
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
[within_subj_stat_ASD] = ft_freqstatistics(cfg, freqFast_mag_ASD{:}, freqSlow_mag_ASD{:})


%plot
cfg = [];
cfg.alpha  = 0.025;
cfg.parameter = 'stat';
cfg.layout = 'neuromag306mag.lay';
ft_clusterplot(cfg, within_subj_stat_NT); 

ft_clusterplot(cfg, within_subj_stat_ASD); 

ft_clusterplot(cfg, within_subj_stat); 

%plot 
saveas(figure(1),[savepath, '/1_results/', 'cluster_based_control.jpeg']);
saveas(figure(2),[savepath, '/1_results/', 'cluster_based_combined.jpeg']);

%% the same for gradiometers
% load freqanalysis data for all groups
for s=1:size(SUBJ_ALL,1)
    
    subj = SUBJ_ALL(s,:); 
    %load data after freqanalysis 
    freq_all = load(strcat(savepath, subj, '/', subj, '_freqanalysis.mat'));
    %after mtmconvol
    freqFast_grad{s}  = freq_all.wvlts_fast_grad;
    freqSlow_grad{s}  = freq_all.wvlts_slow_grad;
end

% load freqanalysis data for control group
for s=1:size(SUBJ_NT,1)
    
    subj = SUBJ_NT(s,:); 
    %load data after freqanalysis 
    freq_all = load(strcat(savepath, subj, '/', subj, '_freqanalysis.mat'));
    %after mtmconvol
    freqFast_grad_NT{s}  = freq_all.wvlts_fast_grad;
    freqSlow_grad_NT{s}  = freq_all.wvlts_slow_grad;
end

% load freqanalysis data for asd group
for s=1:size(SUBJ_ASD,1)
    
    subj = SUBJ_ASD(s,:); 
    %load data after freqanalysis 
    freq_all = load(strcat(savepath, subj, '/', subj, '_freqanalysis.mat'));
    %after mtmconvol
    freqFast_grad_ASD{s}  = freq_all.wvlts_fast_grad;
    freqSlow_grad_ASD{s}  = freq_all.wvlts_slow_grad;
end

%% do cluster-based permutation test
cfg = [];
cfg.channel          = {'MEGGRAD'};
cfg.latency          = [-0.7 0];
cfg.frequency        = [8 12];
cfg.avgoverfreq      = 'yes';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
% specifies with which sensors other sensors can form clusters
cfg_neighb = [];
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, freqFast_grad{1}); %using a distance threshold of 4, there are on average 13.8 neighbours per channel

subj = size (SUBJ_ALL,1);
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

[within_subj_stat_grad] = ft_freqstatistics(cfg, freqFast_grad{:}, freqSlow_grad{:})

%% do cluster-based permutation test for NT
cfg = [];
cfg.channel          = {'MEGGRAD'};
cfg.latency          = [-0.7 0];
cfg.frequency        = [8 12];
cfg.avgoverfreq      = 'yes';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
% specifies with which sensors other sensors can form clusters
cfg_neighb = [];
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, freqFast_grad_NT{1}); %using a distance threshold of 4, there are on average 13.8 neighbours per channel

subj = size (SUBJ_NT,1);
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

[within_subj_stat_grad_NT] = ft_freqstatistics(cfg, freqFast_grad_NT{:}, freqSlow_grad_NT{:})

%% do cluster-based permutation test for ASD
cfg = [];
cfg.channel          = {'MEGGRAD'};
cfg.latency          = [-0.7 0];
cfg.frequency        = [8 12];
cfg.avgoverfreq      = 'yes';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
% specifies with which sensors other sensors can form clusters
cfg_neighb = [];
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, freqFast_grad_ASD{1}); %using a distance threshold of 4, there are on average 13.8 neighbours per channel

subj = size (SUBJ_ASD,1);
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

[within_subj_stat_grad_ASD] = ft_freqstatistics(cfg, freqFast_grad_ASD{:}, freqSlow_grad_ASD{:})


%save stats
filename = strcat(savepath, '1_results/', 'within_subjs_stat_avg_freq.mat');
save(filename, 'within_subj_stat', 'within_subj_stat_NT', 'within_subj_stat_ASD', 'within_subj_stat_grad', 'within_subj_stat_grad_NT', 'within_subj_stat_grad_ASD');

   
%plot
cfg = [];
cfg.alpha  = 0.025;
cfg.parameter = 'stat';
cfg.layout = 'neuromag306planar.lay';%for grad
ft_clusterplot(cfg, within_subj_stat_grad);
ft_clusterplot(cfg, within_subj_stat_grad_NT);
ft_clusterplot(cfg, within_subj_stat_grad_ASD);
%plot 
saveas(figure(2),[savepath, '/1_results/', 'new_cluster_based_grad_combined.jpeg']);
saveas(figure(2),[savepath, '/1_results/', 'cluster_based_grad_control.jpeg']);
saveas(figure(3),[savepath, '/1_results/', 'cluster_based_grad_ASD.jpeg']);
