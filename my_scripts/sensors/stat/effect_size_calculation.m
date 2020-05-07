
close all;

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(path, fieldtripfolder)
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);

realdatapath = '/net/server/data/Archive/aut_gamma/orekhova/KI/SUBJECTS/';
savepath = '/net/server/data/Archive/aut_gamma/orekhova/KI/Scripts_bkp/Shishkina/KI/Results_Alpha_and_Gamma/';

% load subj info  
SUBJ_NT = [ '0101'; '0102'; '0103'; '0104'; '0105'; '0135'; '0136';...  
            '0137'; '0138'; '0140'; '0158'; '0162'; '0163'; '0178';...
            '0179'; '0255'; '0257'; '0348'; '0378'; '0379'; '0384']; 
        
SUBJ_ASD = ['0106'; '0107'; '0139'; '0141'; '0159'; '0160'; '0161';...  
            '0164'; '0253'; '0254'; '0256'; '0273'; '0274'; '0275';...
            '0276'; '0346'; '0347'; '0349'; '0351'; '0358';...
            '0380'; '0381'; '0382'; '0383'];  
        
SUBJ_ALL = [SUBJ_ASD; SUBJ_NT];

freq_fast_isi = cell(1,size(SUBJ_ALL,1));
freq_slow_isi = cell(1,size(SUBJ_ALL,1));

for s=1:size(SUBJ_ALL,1)
    
    subj = SUBJ_ALL(s,:); 
    %load data after freqanalysis 
    load(strcat(savepath, subj, '/', subj, '_alpha_freqanalysis.mat'));
    %after fft
    freq_fast_isi{s} = fft_fast_isi;
    freq_slow_isi{s} = fft_slow_isi;
end

cfg = [];
cfg.keepindividual = 'yes';
grandavg_fast = ft_freqgrandaverage(cfg, freq_fast_isi{:});
grandavg_slow = ft_freqgrandaverage(cfg, freq_slow_isi{:});

% first make the same selection as used in the inferential statistics
cfg = [];
cfg.channel     = 'MEGGRAD';
grandavg_fast_sel = ft_selectdata(cfg, grandavg_fast);
grandavg_slow_sel  = ft_selectdata(cfg, grandavg_slow);

inference = load(strcat(savepath, '1_results/sensor_clusters/', 'cluster_based_stat_grad.mat'));
inference_group = load(strcat(savepath, '1_results/sensor_clusters/', 'cluster_based_between_groups_stat_grad.mat'));
x1 = nan(45,1);
x2 = nan(45,1);

for i=1:45

  % construct a 3-dimensional Boolean array to select the data from this participant
  sel3d = false([45, 204, 40]); 
  sel3d(i,:,:) = inference.stat_all.mask;

  % select the FIC data in the cluster for this participant, represent it as a vector
  tmp = grandavg_fast_sel.powspctrm(sel3d(:));
  % compute the average over the cluster
  x1(i) = mean(tmp);

  % select the FC data in the cluster for this participant, represent it as a vector
  tmp = grandavg_slow_sel.powspctrm(sel3d(:));
  % compute the average over the cluster
  x2(i) = mean(tmp);
end

n1 = length(x1);
n2 = length(x2);

figure; plot([x2 x1]', 'o-'); xlim([0.5 2.5]); xticks([1 2])
xticklabels({'slow','fast'})
legend({'45 subjects'});
title('individual scores, averaged over cluster');

saveas(figure(1),[savepath, '/1_results/sensor_clusters/', 'effect_size_grad.jpeg']);

cohensd = mean(x1-x2) ./ std(x1-x2);
disp(cohensd)