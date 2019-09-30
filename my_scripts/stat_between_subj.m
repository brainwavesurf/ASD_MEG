%the probability distribution of the condition-specific averages is identical for all experimental conditions
clear all;
close all;

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(fieldtripfolder, path);
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);
realdatapath = '/home/a_shishkina/data/KI/SUBJECTS/';
savepath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';
%%
clear all;
close all;
%loading data with powers
SUBJ = ['0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0136'; '0137'; '0138'; '0139'; '0140'; '0141'; '0158'; '0159'; '0160'; '0160'; '0161'; '0162'; '0163'; '0164'; '0178'; '0179'; '0253'; '0254'; '0255'; '0256'; '0257'; '0259'; '0273'; '0274'; '0275'; '0276'; '0277'; '0346'; '0347'; '0348'; '0350'; '0351'; '0357'; '0358'; '0378'; '0380'; '0381'; '0382'; '0383'; '0384']; 

for s = 1: size (SUBJ,1)
    close all
    subj = SUBJ (s,:);
    savepath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';
    allpower = load(strcat(savepath, subj, '/', subj, '_alpha_power.mat'));
    
    pow_slow{s} = allpower.freqPre{1};
    pow_fast{s} = allpower.freqPre{3};
   
end
%% Calculate paired t-test for slow and fast conditions

cfg             = [];
cfg.method      = 'stats'; % using a parametric test
cfg.statistic   = 'paired-ttest'; % using dependent samples
cfg.correctm    = 'fdr'; % no multiple comparisons correction
cfg.alpha       = 0.05;

nsubj           = (size (SUBJ,1));
%cfg.design(1,:) = [1:nsubj 1:nsubj];
cfg.design(1,:) = [ones(1,nsubj)*1 ones(1,nsubj)*2];
%cfg.uvar        = 1; % row of design matrix that contains unit variable (in this case: subjects)
cfg.ivar        = 1; % row of design matrix that contains independent variable (the conditions)

paired_ttest = ft_freqstatistics(cfg, pow_slow{:}, pow_fast{:});


%% Calculate cluster based permutation test
cfg = [];
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
% prepare_neighbours determines what sensors may form clusters
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, pow_slow{1});

nsubj           = (size (SUBJ,1));
cfg.design(1,:) = [1:nsubj 1:nsubj];
cfg.design(2,:) = [ones(1,nsubj)*1 ones(1,nsubj)*2];
cfg.uvar        = 1; % row of design matrix that contains unit variable (in this case: subjects)
cfg.ivar        = 2; % row of design matrix that contains independent variable (the conditions)

[cluster_stat] = ft_freqstatistics(cfg, pow_slow{:}, pow_fast{:});

% calculate the grand average for each condition
cfg = [];
GA_slow       = ft_freqgrandaverage(cfg,pow_slow{:});
GA_fast       = ft_freqgrandaverage(cfg,pow_fast{:});
% "{:}" means to use data from all elements of the variable

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'powspctrm';
GA_fastVSslow = ft_math(cfg,GA_fast,GA_slow);

% In case you have downloaded and loaded the data, ensure stat.cfg.alpha exist
if ~isfield(cluster_stat.cfg,'alpha'); cluster_stat.cfg.alpha = 0.025; end; % stat.cfg.alpha was moved as the downloaded data was processed by an additional FieldTrip function to anonymize the data.

% plot
figure(1);
pos_int = zeros(numel(GA_fastVSslow.label),1);
cfg.highlight = 'on';
cfg.highlightchannel = find(pos_int);
cfg.commentpos = 'title';
cfg.layout = 'neuromag306planar_helmet.mat';
ft_topoplotER(cfg, GA_fastVSslow);

saveas(figure(1),[savepath, '1_results_plot/', 'cluster_based_PSD.jpeg']);

filename = strcat(savepath, '1_results_plot/', 'paired_ttest.mat');
save (filename, 'paired_ttest', 'cluster_stat');