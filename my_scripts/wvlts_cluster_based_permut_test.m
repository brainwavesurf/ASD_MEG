%% Cluster-based permutation test for within-subjects experiments after wvlts freqanalysis
%(http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/#within-subjects-experiments)
%for magnetometers and gradiometers
%stats for 5to30 Hz, [-0.8 0] interstimuli interval
%clusterplots for frequencies with significant fast vs slow sensory stimuli intensity
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

% Define sensors here: MAG or GRAD
channels     = 'megmag'; %'megplanar';% 
layout_mat      = 'neuromag306mag_helmet.mat'; %'neuromag306planar_helmet.mat';%
layout_lay      = 'neuromag306mag.lay'; %'neuromag306planar.lay'; %
template     = 'neuromag306mag_neighb.mat';  %'neuromag306planar_neighb.mat'; %
screensize = get( groot, 'Screensize' );

Nsub = size(SUBJ_ALL, 1);
%% load freqanalysis data for all groups
for s = 1:Nsub
    
    subj = SUBJ_ALL(s,:); 
    %load data after freqanalysis 
    freq_all = load(strcat(savepath, subj, '/', subj, '_freqanalysis.mat'));;
    %after mtmconvol
    freqFast_grad{s}  = freq_all.wvlts_fast_mag; %freq_all.wvlts_fast_grad; %
    freqSlow_grad{s}  = freq_all.wvlts_slow_mag; %freq_all.wvlts_slow_grad; %
end

cfg = [];
cfg.method      = 'template'; %'distance'; % try as well
cfg.template    = template;               % specify type of template
cfg.layout      = layout_mat;               % specify layout of sensors*
cfg.feedback    = 'yes';                             % show a neighbour plot
neighbours      = ft_prepare_neighbours(cfg, freqFast_mag{1}); 
%neighbours      = ft_prepare_neighbours(cfg, freqFast_grad{1}); 
%%

cfg=[];
cfg.channel     = channels; % 'MEG'; %'megplanar'
cfg.neighbours  = neighbours;
cfg.parameter   = 'powspctrm'; 
cfg.frequency   = [5 30];
cfg.latency     = [-0.8 0];
cfg.avgovertime = 'yes';
cfg.method      = 'montecarlo'; % 'stats';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.minnbchan = 1;

cfg.correctm    = 'cluster'; %'fdr'
cfg.correcttail = 0.025;
cfg.numrandomization = 10000; 
cfg.clusterstatistic = 'maxsum'; % 'maxsum', 'maxsize', 'wcm'-'weighted cluster mass'; (default = 'maxsum')
cfg.clusterthreshold = 'parametric'; % method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric')
cfg.clusteralpha     = 0.05; % for either parametric or nonparametric thresholding per tail (default = 0.05)
cfg.clustertail      = 0; %-1, 1 or 0 (default = 0)


cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable %The “Independent variable” codes the condition numbers. 
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number       %The “Unit of observation variable” corresponds to the subject number 

[within_subj_stat_mag] = ft_freqstatistics(cfg, freqFast_mag{:}, freqSlow_mag{:});
%[within_subj_stat_grad] = ft_freqstatistics(cfg, freqFast_grad{:}, freqSlow_grad{:});

% %save stats
filename = strcat(savepath, '1_results/', 'within_subjs_stat_wvlts_mag.mat');
save(filename, 'within_subj_stat_mag');
%filename = strcat(savepath, '1_results/', 'within_subjs_stat_wvlts_grad.mat');
%save(filename, 'within_subj_stat_grad');
%% clusterplot
cfg = [];
cfg.zlim = [-4 4]; % T-values
cfg.alpha = 0.05;
cfg.parameter = 'stat';
cfg.channel = channels; 
cfg.layout = layout_lay;
cfg.highlightcolorpos =  [1 0 0]; %color of highlight marker for positive clusters (default = [0 0 0])
ft_clusterplot(cfg, within_subj_stat_mag);
%ft_clusterplot(cfg, within_subj_stat_grad);
colorbar

 for i = 1:2
       
   saveas(figure(i),[savepath, '1_results/', num2str(i), '_wvlts_freq_analysis_mag_round.jpeg']);
   %saveas(figure(i),[savepath, '1_results/', num2str(i), '_wvlts_freq_analysis_grad.jpeg']);
   end

%% plot clusters, join frequencies
% NB: inspect number of significant clusters and correct subplot size accordingly

cfg =[];
cfg.layout = layout_lay;
layout = ft_prepare_layout(cfg, []);
pos_fig2 = [10 10 screensize(3)/4*3 screensize(4)/2];    
hh = figure('Position',pos_fig2);

%NN=number_of_significant_clusters
NN=1;
for n=1:NN
    clusterN = n;
    [ch,fr] = find(within_subj_stat_mag.posclusterslabelmat==clusterN);
    Fmin = within_subj_stat_mag.freq(min(fr));
    Fmax = within_subj_stat_mag.freq(max(fr));
    CHAN = unique(ch);
   
    subplot(1,NN,n)
    ft_plot_layout (layout, 'label', 'yes', 'chanindx',  CHAN );
    title ({['size:', num2str(size(fr)),] ['Freq:', num2str(Fmin), '-', num2str(Fmax), 'Hz'] ['sig: p=', num2str(within_subj_stat.posclusters(n).prob)]}); 
end

suptitle ('clusters for All subj, Magnetometers, 5 to 30 Hz, [-0.8 0] sec')
%suptitle ('clusters for All subj, Gradiometers, 5 to 30 Hz, [-0.8 0] sec')

saveas(figure(1),[savepath, '/1_results/', 'mag_clusters_location.jpeg']);
%saveas(figure(1),[savepath, '/1_results/', 'grad_clusters_location.jpeg']);