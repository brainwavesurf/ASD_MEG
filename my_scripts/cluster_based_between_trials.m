%% Cluster-based permutation test for between-trial experiments 
%(http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/#within-subjects-experiments)
% my questions are commented and boundared by "===", also I added info about stuctures 
% and as a result I got:
%Error using ft_clusterplot (line 208)
%no clusters present with a p-value lower than the specified alpha, nothing to plot
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
SUBJ = ['0102'];

for s=1: size (SUBJ,1)
    close all
    subj = SUBJ (s,:); 
    savemegto = strcat(savepath, subj);
    epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
    
    %% load group preceding events in order to select the epochs later on
    load ([ savemegto, '/', subj, '_info.mat'])
    
    slow_ind = find(allinfo.prev_stim_type==2);
    fast_ind = find(allinfo.prev_stim_type==8);
 
    %% Load unfiltered epochs and divide them according to preceding conditions
    ep_fiff_file = strcat(epofolder, subj, '-noerror-lagcorrected-epo.fif')
    hdr = ft_read_header(ep_fiff_file);
    
    cfg           = [];  
    cfg.dataset   = ep_fiff_file;
    cfg.channel   = {'MEG'};
    epochs        = ft_preprocessing(cfg);
    
    cfg           = [];
    cfg.trials    = slow_ind;
    slow_epochs   = ft_selectdata(cfg, epochs); %extraxt trials after slow stimuli
    
    cfg.trials    = fast_ind;
    fast_epochs   = ft_selectdata(cfg, epochs); %extraxt trials after fast stimuli
    
    %% Calculation of the planar gradient and time-frequency analysis
    cfg = [];
    cfg.planarmethod     = 'sincos';
    % prepare_neighbours determines with what sensors the planar gradient is computed
    cfg_neighb.method    = 'distance';
    cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, fast_epochs);
    %Warning: Please be aware of the different sensor types in neuromag306 system, see
    %http://www.fieldtriptoolbox.org/faq/why_are_there_multiple_neighbour_templates_for_the_neuromag306_system 
    %======================================================================
    % ??? How could I change my data to avoid this warning and how important it is
    %======================================================================
    
    dataF_planar  = ft_megplanar(cfg, slow_epochs);
    dataS_planar  = ft_megplanar(cfg, fast_epochs);
    
    %Error using ft_checkdata (line 574)
    %This function requires data with an 'ctf151', 'ctf275', 'bti148', 'bti248', 'itab153', 'yokogawa160' or 'yokogawa64' sensor array.
    %Error in ft_megplanar (line 97)
    %data = ft_checkdata(data, 'datatype', {'raw' 'freq'}, 'feedback', 'yes', 'hassampleinfo', 'yes', 'ismeg', 'yes', 'senstype',
    %{'ctf151', 'ctf275', 'bti148', 'bti248', 'itab153', 'yokogawa160', 'yokogawa64'});
    %======================================================================
    % ???What method I should use for Neuromag Elekta MEG system which is
    % not presented here (I skip this step for further analysis and do 
    %ft_combineplanar on preprocessed data instead of data after ft_megplanar. How to solve this problem? 
    %======================================================================
    
    %% Calculate the TFRs for the two experimental conditions (Fast and Slow) for each subject

    cfg = [];
    cfg.output     = 'pow';
    cfg.channel    = 'MEG';
    cfg.method     = 'mtmconvol';
    cfg.taper      = 'hanning';
    cfg.foi        = [5:20];
    cfg.toi        = [-1:0.05:1.2];
    cfg.t_ftimwin  = 7./cfg.foi; %7 cycles
    cfg.keeptrials = 'yes';

    freqFast = ft_freqanalysis(cfg, fast_epochs);
    freqSlow = ft_freqanalysis(cfg, slow_epochs);
    
    
    %% calculate the combined planar gradient and copy the gradiometer structure in the new datasets
    cfg = [];
    cfg.latency = [-0.8 0];
    freqFast = ft_selectdata(cfg, freqFast);
    freqSlow = ft_selectdata(cfg, freqSlow);
    
    cfg = [];
    freqSlow_planar = ft_combineplanar(cfg, freqSlow);
    freqFast_planar = ft_combineplanar(cfg, freqFast);

    freqSlow_planar.grad = slow_epochs.grad;
    freqFast_planar.grad = fast_epochs.grad;
    
    %freqFast_planar = 
    %struct with fields:
    %    label: {204×1 cell}
    %   dimord: 'rpt_chan_freq_time'
    %     freq: [1×16 double]
    %     time: [1×17 double]
    %powspctrm: [63×204×16×17 double]
    %cumtapcnt: [63×16 double]
    %     grad: [1×1 struct]
    %      cfg: [1×1 struct]
    
    %% Do between trials stats
    cfg = [];
    cfg.channel          = {'MEG'};
    cfg.latency          = [-0.8 0];
    cfg.frequency        = [5 20];
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'ft_statfun_indepsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan        = 2;
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.alpha            = 0.025;
    cfg.numrandomization = 500;
    
    cfg_neighb.method    = 'distance';
    cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, fast_epochs);

    design = zeros(1,size(freqSlow_planar.powspctrm,1) + size(freqFast_planar.powspctrm,1));
    design(1,1:size(freqSlow_planar.powspctrm,1)) = 1;
    design(1,(size(freqSlow_planar.powspctrm,1)+1):(size(freqSlow_planar.powspctrm,1)+...
    size(freqFast_planar.powspctrm,1))) = 2;
    
    cfg.design           = design;
    cfg.ivar             = 1;
    
    [stat] = ft_freqstatistics(cfg, freqSlow_planar, freqFast_planar);
    
%% stat results

% stat = 
% 
%   struct with fields:
% 
%                    prob: [204×16×17 double]
%             posclusters: [1×10 struct]
%     posclusterslabelmat: [204×16×17 double]
%         posdistribution: [1×500 double]
%             negclusters: [1×7 struct]
%     negclusterslabelmat: [204×16×17 double]
%         negdistribution: [1×500 double]
%                 cirange: [204×16×17 double]
%                    mask: [204×16×17 logical]
%                    stat: [204×16×17 double]
%                     ref: [204×16×17 double]
%                  dimord: 'chan_freq_time'
%                    freq: [1×16 double]
%                    grad: [1×1 struct]
%                   label: {204×1 cell}
%                    time: [1×17 double]
%                     cfg: [1×1 struct]


    %% do average over trials
    cfg = [];
    freqFast_avg = ft_freqdescriptives(cfg, freqFast_planar);
    freqSlow_avg = ft_freqdescriptives(cfg, freqSlow_planar);
    
%     freqSlow_avg = 
%   struct with fields:
% 
%        dimord: 'chan_freq_time'
%          freq: [1×16 double]
%         label: {204×1 cell}
%     powspctrm: [204×16×17 double]
%          time: [1×17 double]
%          grad: [1×1 struct]
%     cumtapcnt: [68×16 double]
%           cfg: [1×1 struct]
    %% add the raw effect (Fast-Slow) to the obtained stat structure and plot the largest cluster overlayed on the raw effect.
    stat.raweffect = freqFast_avg.powspctrm - freqSlow_avg.powspctrm; %raweffect: [204×16×17 double]
    
    cfg = [];
    cfg.alpha  = 0.025;
    cfg.parameter = 'raweffect';
    cfg.zlim   = [-1e-27 1e-27];
    cfg.layout = 'neuromag306mag.lay';
    ft_clusterplot(cfg, stat);
    
    %%Error using ft_clusterplot (line 159)
    %unsupported dimord chan_unknown_time
    %======================================================================
    %to solve this problem a do averaging for frequencies. Is it correctly?
    %I used advice from https://mailman.science.ru.nl/pipermail/fieldtrip/2016-May/010415.html
    %Are there other ways to fix it?
    %======================================================================
    raweffect = freqFast_avg; % for apropriate structure
    raweffect.powspctrm = freqFast_avg.powspctrm - freqSlow_avg.powspctrm; %raweffect: [204×1×17 double]
    
    cfg  = [];
    cfg.avgoverfreq = 'yes';
    raweffect = ft_selectdata(cfg, raweffect);
    stat.raweffect = raweffect.powspctrm;
    
    cfg  = [];
    cfg.avgoverfreq = 'yes';
    stat = ft_selectdata(cfg, stat);
    
%                    prob: [204×1×17 double]
%             posclusters: [1×10 struct]
%     posclusterslabelmat: [204×1×17 double]
%         posdistribution: [1×500 double]
%             negclusters: [1×7 struct]
%     negclusterslabelmat: [204×1×17 double]
%         negdistribution: [1×500 double]
%                 cirange: [204×1×17 double]
%                    mask: [204×1×17 double]
%                    stat: [204×1×17 double]
%                     ref: [204×1×17 double]
%                    freq: 12.4896
%                    grad: [1×1 struct]
%                   label: {204×1 cell}
%                    time: [1×17 double]
%                     cfg: [1×1 struct]
%               raweffect: [204×1×17 double]
%                  dimord: 'chan_freq_time'

    cfg = [];
    cfg.alpha  = 0.025;
    cfg.parameter = 'raweffect';
    %cfg.zlim   = [-1e-27 1e-27];
    cfg.layout = 'neuromag306mag.lay';
    ft_clusterplot(cfg, stat);

    %Error using ft_clusterplot (line 208)
    %no clusters present with a p-value lower than the specified alpha, nothing to plot
end 


