%the probability distribution of the condition-specific averages is identical for all experimental conditions
%%
clear all;
close all;
%loading data with powers
SUBJ = ['0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0136'; '0137'; '0138'; '0139'; '0140'];
% posterior sensors, you may use them later for extracting alpha peak power
post_ch = {'MEG1932',  'MEG1922', 'MEG2042',  'MEG2032',  'MEG2112', 'MEG2122',  'MEG2342', 'MEG2332',  'MEG1732', 'MEG1942', 'MEG1912', 'MEG2012', 'MEG2022', 'MEG2312', 'MEG2322', 'MEG2512',...
%           'MEG1933',  'MEG1923', 'MEG2043',  'MEG2033',  'MEG2113', 'MEG2123',  'MEG2343', 'MEG2333',  'MEG1733', 'MEG1943', 'MEG1913', 'MEG2013', 'MEG2023', 'MEG2313', 'MEG2323', 'MEG2513'};

%% grand averaging
cfg = [],
tfr_avg_ = ft_freqgrandaverage([], tfr_slow);


cfg = [];
tfr_app = ft_appendfreq(cfg, tfr_slow, tfr_fast);

for s = 1: size (SUBJ,1)
    close all
    subj = SUBJ (s,:);
    savepath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';
    allpower = load(strcat(savepath, subj, '/', subj, '_alpha_power.mat'));
    
    
    %allsubj_slow{s} = allpower.powPre{1};
    %allsubj_fast{s} = allpower.powPre{3};
    
    %[h,p] = ttest(allsubj_slow{s}, allsubj_fast{s});
    %p_val{s} = p; 
    
   
end


tfr_slow(s) = allpower.freqPre{1};
tfr_fast(s) = allpower.freqPre{3};


 cfg             = [];
cfg.method      = 'stats'; % using a parametric test
cfg.statistic   = 'depsamplesT'; % using dependent samples
cfg.correctm    = 'fdr'; % no multiple comparisons correction
cfg.alpha       = 0.05;

nsubj           = numel(size (SUBJ,1));
cfg.design(1,:) = [1:nsubj 1:nsubj];
cfg.design(2,:) = [ones(1,nsubj)*1 ones(1,nsubj)*2];
cfg.uvar        = 1; % row of design matrix that contains unit variable (in this case: subjects)
cfg.ivar        = 2; % row of design matrix that contains independent variable (the conditions)

stat{s} = ft_freqstatistics(cfg, allpower.freqPre{1}, allpower.freqPre{3});

