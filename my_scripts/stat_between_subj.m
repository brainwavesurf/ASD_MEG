%%
clear all;
close all;
%loading data with powers
SUBJ = ['0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0136'; '0137'; '0138'; '0139'; '0140'];

for s = 1: size (SUBJ,1)
    close all
    subj = SUBJ (s,:);
    savepath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';
    allpower = load(strcat(savepath, subj, '/', subj, '_alpha_power.mat'));
    
    allsubj_slow{s} = allpower.powPre{1};
    allsubj_fast{s} = allpower.powPre{3};
    
    [h,p] = ttest(allsubj_slow{s}, allsubj_fast{s});
    p_val{s} = p; 
end

fdr = mafdr(p); 

cfg             = [];
cfg.method      = 'stats'; % using a parametric test
cfg.statistic   = 'depsamplesT'; % using dependent samples
%cfg.correctm    = 'fdr'; % no multiple comparisons correction
cfg.alpha       = 0.05;

nsubj           = numel(size (SUBJ,1));
cfg.design(1,:) = [1:nsubj 1:nsubj];
cfg.design(2,:) = [ones(1,nsubj)*1 ones(1,nsubj)*2];
cfg.uvar        = 1; % row of design matrix that contains unit variable (in this case: subjects)
cfg.ivar        = 2; % row of design matrix that contains independent variable (the conditions)

stat = ft_freqstatistics(cfg, allsubj_slow{:}, allsubj_fast{:});
[h,p] = ttest(allsubj_slow{1}, allsubj_fast{1})
p_val{s} = p 

