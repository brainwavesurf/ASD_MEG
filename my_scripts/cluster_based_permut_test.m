%% Cluster-based permutation test for within-subjects experiments 
%(http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/#within-subjects-experiments)
% Reading the data for trials in slow and fast conditions
% Do time-frequency transformation with wavelets for [-1 1.2] time interval
% Do the grand averages of the TFRs for the fast and slow conditions
% Do group-level stats for interval of interest [-0.8 0] and frequency range [5 20]
%==========================================================================
% error!!!   
%Error using ft_statfun_depsamplesT (line 84)
%The data must contain at least two units (usually subjects).
%Error using ft_statistics_montecarlo (line 244)
%could not determine the parametric critical value for clustering
%Error in ft_freqstatistics (line 193)
%[stat, cfg] = statmethod(cfg, dat, design);
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
SUBJ = ['0102';'0103'; '0104'];

for s=1: size (SUBJ,1)
    close all
    subj = SUBJ (s,:); 
    savemegto = strcat(savepath, subj);
    epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
    
    %% load group preceding events in order to select the epochs later on
    load ([ savemegto, '/', subj, '_info.mat'])
    
    slow_ind = find(allinfo.prev_stim_type==2);
    fast_ind = find(allinfo.prev_stim_type==8);
 
    %% Load unfiltered epochs 
    ep_fiff_file = strcat(epofolder, subj, '-noerror-lagcorrected-epo.fif')
    hdr = ft_read_header(ep_fiff_file);
    
    cfg           = [];  
    cfg.dataset   = ep_fiff_file;
    cfg.channel   = {'MEG'};
    epochs        = ft_preprocessing(cfg);
    
    cfg           = [];
    cfg.trials    = slow_ind;
    fast_epochs   = ft_selectdata(cfg, epochs);
    
    cfg.trials    = fast_ind;
    slow_epochs   = ft_selectdata(cfg, epochs);
    
    %% Calculate the TFRs for the two experimental conditions (Fast and Slow) for each subject
    cfg = [];
    cfg.channel    = 'MEG';
    cfg.method     = 'wavelet';
    cfg.width      = 7;
    cfg.output     = 'pow';
    cfg.foi        = 5:20;
    cfg.toi        = -1:0.05:1.2;

    freqFast{s}  = ft_freqanalysis(cfg, fast_epochs);
    freqSlow{s}  = ft_freqanalysis(cfg, slow_epochs);
    
    
    %Plot the results
    %cfg = [];
    %cfg.showlabels   = 'yes';
    %cfg.layout       = 'neuromag306planar_helmet.mat';
    %figure
    %ft_multiplotTFR(cfg, freqFast)
end
%% Calculate the grand averages of the TFRs for the fast and slow conditions

cfg = [];
cfg.foilim = [5 20];
cfg.toilim = [-0.8 0];
cfg.channel = 'meg';
slow_avg = ft_freqgrandaverage(cfg, freqSlow{:});
fast_avg = ft_freqgrandaverage(cfg, freqFast{:});

%Averaging of grads info because of warning "discarding gradiometer information because it cannot be averaged" 
AvgGrad_fast = ft_average_sens([freqFast{1}.grad, freqFast{2}.grad, freqFast{3}.grad]);
fast_avg.grad = AvgGrad_fast;

AvgGrad_slow = ft_average_sens([freqSlow{1}.grad, freqSlow{2}.grad, freqSlow{3}.grad]);
slow_avg.grad = AvgGrad_slow;

%% do cluster-based permutation test
cfg = [];
cfg.channel          = {'MEG'};
cfg.latency          = [-0.8 0];
cfg.frequency        = [5 20];
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
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, fast_avg); %using a distance threshold of 4, there are on average 13.8 neighbours per channel

subj = size (SUBJ,1);
design = zeros(2,2*subj);
design(1,1:subj)        = 1;
design(1,subj+1:2*subj) = 2;
for i = 1:subj
design(2,i) = i;
end
for i = 1:subj
design(2,subj+i) = i;
end

cfg.design   = design;

cfg.uvar     = 1;
cfg.ivar     = 2;

[stat] = ft_freqstatistics(cfg, fast_avg, slow_avg)

