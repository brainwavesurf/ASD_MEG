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

filename = strcat(savepath, '1_results_plot/', 'paired_ttest.mat');
save (filename, 'paired_ttest');