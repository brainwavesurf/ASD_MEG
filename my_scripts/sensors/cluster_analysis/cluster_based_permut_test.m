%% Cluster-based permutation test for within-subjects experiments 
%(http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/#within-subjects-experiments)
%for magnetometers
%for gradiometers
%all subjects, NT, ASD groups
%%
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
%% load freqanalysis data for all groups

freq_fast_isi = cell(1,size(SUBJ_ALL,1));
freq_slow_isi = cell(1,size(SUBJ_ALL,1));
nsub_all = size(SUBJ_ALL,1);

for s=1:size(SUBJ_ALL,1)
    
    subj = SUBJ_ALL(s,:); 
    %load data after freqanalysis 
    load(strcat(savepath, subj, '/', subj, '_alpha_freqanalysis.mat'));
    %after fft
    freq_fast_isi{s} = fft_fast_isi;
    freq_slow_isi{s} = fft_slow_isi;
end


%% load freqanalysis data for control group
freq_fast_isi_NT = cell(1,size(SUBJ_NT,1));
freq_slow_isi_NT = cell(1,size(SUBJ_NT,1));
nsub_NT = size(SUBJ_NT,1);
for s=1:size(SUBJ_NT,1)
    
    subj = SUBJ_NT(s,:); 
    %load data after freqanalysis 
    load(strcat(savepath, subj, '/', subj, '_alpha_freqanalysis.mat'));
    %after fft
    freq_fast_isi_NT{s} = fft_fast_isi;
    freq_slow_isi_NT{s} = fft_slow_isi;
end

%% load freqanalysis data for asd group
freq_fast_isi_ASD = cell(1,size(SUBJ_ASD,1));
freq_slow_isi_ASD = cell(1,size(SUBJ_ASD,1));
nsub_ASD = size(SUBJ_ASD,1);
for s=1:size(SUBJ_ASD,1)
    
    subj = SUBJ_ASD(s,:); 
    %load data after freqanalysis 
    load(strcat(savepath, subj, '/', subj, '_alpha_freqanalysis.mat'));
    %after fft
    freq_fast_isi_ASD{s} = fft_fast_isi;
    freq_slow_isi_ASD{s} = fft_slow_isi;
end

%mag 
channels     = 'megmag'; 
layout1      = 'neuromag306mag_helmet.mat'; 
layout2      = 'neuromag306mag.lay'; 
template     = 'neuromag306mag_neighb.mat';  

% %grad
% channels     = 'megplanar'; 
% layout1      = 'neuromag306planar_helmet.mat'; 
% layout2      = 'neuromag306planar.lay'; 
% template     = 'neuromag306planar_neighb.mat'; 
%% do cluster-based permutation test for ALL
cfg = [];
cfg.method      = 'template'; 
cfg.template    = template;              
cfg.layout      = layout2;             
cfg.feedback    = 'yes';                             
neighbours      = ft_prepare_neighbours(cfg, freq_fast_isi{1});
neighbours_NT   = ft_prepare_neighbours(cfg, freq_fast_isi_NT{1});
neighbours_ASD  = ft_prepare_neighbours(cfg, freq_fast_isi_ASD{1});

cfg=[];
cfg.channel     = channels; 
cfg.parameter   = 'powspctrm';
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.minnbchan   = 1;
cfg.correctm    = 'cluster'; %'fdr'
cfg.correcttail = 0.025;
cfg.numrandomization = 10000; 
cfg.clusterstatistic = 'maxsum'; 
cfg.clusterthreshold = 'parametric'; 
cfg.clusteralpha     = 0.05; 
cfg.clustertail      = 0; 

cfg.design(1,1:2*nsub_all)  = [ones(1,nsub_all) 2*ones(1,nsub_all)];
cfg.design(2,1:2*nsub_all)  = [1:nsub_all 1:nsub_all];
cfg.ivar                = 1; 
cfg.uvar                = 2; 
cfg.neighbours  = neighbours;
[stat_all] = ft_freqstatistics(cfg, freq_fast_isi{:}, freq_slow_isi{:});

cfg.design = [];
cfg.design(1,1:2*nsub_NT)  = [ones(1,nsub_NT) 2*ones(1,nsub_NT)];
cfg.design(2,1:2*nsub_NT)  = [1:nsub_NT 1:nsub_NT];
cfg.ivar                = 1; 
cfg.uvar                = 2; 
cfg.neighbours  = neighbours_NT;
[stat_NT] = ft_freqstatistics(cfg, freq_fast_isi_NT{:}, freq_slow_isi_NT{:});

cfg.design = [];
cfg.design(1,1:2*nsub_ASD)  = [ones(1,nsub_ASD) 2*ones(1,nsub_ASD)];
cfg.design(2,1:2*nsub_ASD)  = [1:nsub_ASD 1:nsub_ASD];
cfg.ivar                = 1; 
cfg.uvar                = 2; 
cfg.neighbours  = neighbours_ASD;
[stat_ASD] = ft_freqstatistics(cfg, freq_fast_isi_ASD{:}, freq_slow_isi_ASD{:});

%save stats
filename = strcat(savepath, '1_results/sensor_clusters/', 'cluster_based_stat_mag.mat');
save(filename, 'stat_all', 'stat_NT', 'stat_ASD')

% filename = strcat(savepath, '1_results/sensor_clusters/', 'cluster_based_stat_grad.mat');
% save(filename, 'stat_all', 'stat_NT', 'stat_ASD')



%plot
cfg = [];
cfg.alpha  = 0.025;
cfg.parameter = 'stat';
cfg.layout = layout2;
ft_clusterplot(cfg, stat_all); 
ft_clusterplot(cfg, stat_NT); 
ft_clusterplot(cfg, stat_ASD); 
 
saveas(figure(1),[savepath, '/1_results/sensor_clusters/', 'cluster_all_mag.jpeg']);
saveas(figure(2),[savepath, '/1_results/sensor_clusters/', 'cluster_NT_mag.jpeg']);
saveas(figure(1),[savepath, '/1_results/sensor_clusters/', 'cluster_ASD_mag.jpeg']);

saveas(figure(2),[savepath, '/1_results/sensor_clusters/', 'cluster_all_grad_2.jpeg']);
saveas(figure(1),[savepath, '/1_results/sensor_clusters/', 'cluster_NT_grad.jpeg']);
saveas(figure(1),[savepath, '/1_results/sensor_clusters/', 'cluster_ASD_grad.jpeg']);


load(strcat(savepath, '1_results/sensor_clusters/', 'cluster_based_stat_mag.mat'));
%load(strcat(savepath, '1_results/sensor_clusters/', 'cluster_based_stat_grad.mat'));
%NN=number_of_significant_clusters
NN=2;
cfg=[];cfg.layout=layout2;
layout = ft_prepare_layout(cfg, []);
for n=1:NN
    clusterN=n;
    [ch, fr] = find(stat_all.posclusterslabelmat == clusterN);
    F1 = stat_all.freq(min(fr));
    F2 = stat_all.freq(max(fr));
    CH = unique(ch);
   
    subplot(2,1,n)
    ft_plot_layout(layout, 'label', 'yes', 'chanindx',  CH);
    title (['size:', num2str(size(fr)),', Freq:', num2str(F1), '-', num2str(F2), 'Hz, sig: p=', num2str(stat_all.posclusters(n).prob)]); 
   %%
end

suptitle ('clusters for all subj, Magnetometers, 2 to 40 Hz')
%suptitle ('clusters for ASD subj, Gradiometers, 2 to 40 Hz')

saveas(figure(1),[savepath, '/1_results/sensor_clusters/', 'cluster_channels_mag_all.jpeg']);
%saveas(figure(1),[savepath, '/1_results/sensor_clusters/', 'cluster_channels_grad_ASD.jpeg']);