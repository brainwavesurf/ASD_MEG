%%
clear
close all
screensize = get( groot, 'Screensize' );

fieldtripfolder = '/mnt/home/a_shishkina/fieldtrip/';

path(path, fieldtripfolder)
ft_defaults;

path( strcat (fieldtripfolder, 'external/mne/'), path);
path(path,'/mnt/home/a_shishkina/externals/pptx/');
path(path,'/mnt/home/a_shishkina/projects/asd_meg/0101/MEGdata/FT/'); %for fdr_bh.m

realdatapath = '/mnt/home/a_shishkina/data/KI/SUBJECTS/'; 
DataPath = '/mnt/home/a_shishkina/data/KI/FT_beamformer/';

tapsmofrq = 2; 
gridres = 6; 
lambda='5%';

alpharange = [5 20];
alphamax = [7 14];

Nvoxels = 25;

filt = 'yes';
lpfreq=alpharange(2)+5;
hpfreq=[];

megfolder = strcat( '/meg', num2str(gridres), 'mm_linwarp_tapsmofrq', num2str(tapsmofrq), '_lcmv_alpha_', '/');

SUBJ = ['0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0136'; '0137'; '0138'; '0139'; '0140'; '0141'; '0158'; '0159'; '0160'; '0161'; '0162'; '0163'; '0164'; '0178'; '0179'; '0253'; '0254'; '0255'; '0256'; '0257'; '0259'; '0273'; '0274'; '0275'; '0276'; '0277'; '0346'; '0347'; '0348'; '0350'; '0351'; '0357'; '0358'; '0378'; '0380'; '0381'; '0382'; '0383'; '0384'; '0385']; % 

subj = SUBJ(1,:);

%% Load clean events
load (strcat (realdatapath, subj, '/ICA_nonotch_crop/', subj, '_clean_events.mat'));

%% Load raw data
fiff_file = strcat(realdatapath, subj, '/ICA_nonotch_crop/', subj, '_rings_ICA_raw.fif');
hdrraw = ft_read_header(fiff_file);

cfg =[];
layout = ft_prepare_layout(cfg, fiff_file);

pre = -0.8* hdrraw.Fs;
post = 1.2* hdrraw.Fs;

first = round(cast(hdrraw.orig.raw.first_samp, 'double'));
events(:,1) =  events(:,1)-first;

trl = [];
for i=1:size (events,1)
     trl(i, 1)=(events(i,1)+pre) ; 
     trl(i, 2)=(events(i,1)+post) ; 
     trl(i, 3)= -1.0*hdrraw.Fs ; % offset
     trl(i, 4) = events(i,3); % stimulus_value;
end

%% extract data and epochs from the raw
cfg = [];
cfg.dataset = fiff_file;
cfg.trl = trl;
cfg = ft_definetrial(cfg);

cfg.channel = 'meg';
cfg.dftfilter = 'yes';
cfg.dftfreq = [50 100];
cfg.demean = 'yes';
cfg.lpfilter = 'yes';
cfg.lpfreq = lpfreq;
epochs = ft_preprocessing(cfg);

%%  select epochs according to events

evslow = find(events(:,3)==2);
evfast = find(events(:,3)==8);

%% select epochs slow, medium, fast
cfg = [];
cfg.trials =  evslow; 
eposlow = ft_selectdata(cfg, epochs);

cfg = [];
cfg.trials =  evfast; 
epofast = ft_selectdata(cfg, epochs);

cfg = [];
cfg.channel=eposlow.label;
cfg.latency = [-0.8 0];
cfg.keeptrials = 'yes';
avgslow = ft_timelockanalysis(cfg, eposlow);

cfg = [];
cfg.channel=epofast.label;
cfg.latency = [-0.8 0];
cfg.keeptrials = 'yes';
avgfast = ft_timelockanalysis(cfg, epofast);

%% cluster based
cfg = [];
cfg.method = 'montecarlo';       % use the Monte Carlo Method to calculate the significance probability
cfg.statistic = 'ft_statfun_indepsamplesT'; % use the independent samples T-statistic as a measure to
                               % evaluate the effect at the sample level
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;         % alpha level of the sample-specific test statistic that
                               % will be used for thresholding
cfg.clusterstatistic = 'maxsum'; % test statistic that will be evaluated under the
                               % permutation distribution.
cfg.minnbchan = 2;               % minimum number of neighborhood channels that is
                               % required for a selected sample to be included
                               % in the clustering algorithm (default=0).

cfg.tail = 0;                    % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail = 0;
cfg.alpha = 0.025;               % alpha level of the permutation test
cfg.numrandomization = 100;      % number of draws from the permutation distribution

design = zeros(1,size(avgslow.trial,1) + size(avgfast.trial,1));
design(1,1:size(avgslow.trial,1)) = 1;
design(1,(size(avgslow.trial,1)+1):(size(avgslow.trial,1) + size(avgfast.trial,1)))= 2;

cfg.design = design;             % design matrix
cfg.ivar  = 1;

%% neighbours
cfg_neighb        = [];
cfg_neighb.method = 'distance';
neighbours        = ft_prepare_neighbours(cfg_neighb, epofast);

cfg.neighbours    = neighbours;  % the neighbours specify for each sensor with
                                 % which other sensors it can form clusters
cfg.channel       = {'MEG'};     % cell-array with selected channel labels
cfg.latency       = [-0.8 0];  

%% stat
[stat] = ft_timelockstatistics(cfg, avgslow, avgfast);

%% plot
cfg = [];
cfg.channel=eposlow.label;
cfg.latency = [-0.8 0];
avgslow = ft_timelockanalysis(cfg, eposlow);

cfg = [];
cfg.channel=epofast.label;
cfg.latency = [-0.8 0];
avgfast = ft_timelockanalysis(cfg, epofast);

cfg  = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
raweffectSLOWvsFAST = ft_math(cfg,avgslow,avgfast);

%%
pos_cluster_pvals = [stat.posclusters(:).prob];

% Then, find which clusters are significant, outputting their indices as held in stat.posclusters
% In case you have downloaded and loaded the data, ensure stat.cfg.alpha exist
if ~isfield(stat.cfg,'alpha'); stat.cfg.alpha = 0.025; end; % stat.cfg.alpha was moved as the downloaded data was processed by an additional FieldTrip function to anonymize the data.

pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
% (stat.cfg.alpha is the alpha level we specified earlier for cluster comparisons; In this case, 0.025)
% make a boolean matrix of which (channel,time)-pairs are part of a significant cluster
pos = ismember(stat.posclusterslabelmat, pos_signif_clust);

% and now for the negative clusters...
neg_cluster_pvals = [stat.negclusters(:).prob];
neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
neg = ismember(stat.negclusterslabelmat, neg_signif_clust);

%% manipulate clasters
pos = stat.posclusterslabelmat == 1; % or == 2, or 3, etc.
neg = stat.negclusterslabelmat == 1;

%%
timestep = 0.05; % timestep between time windows for each subplot (in seconds)
sampling_rate = epofast.fsample; % Data has a temporal resolution of 300 Hz
sample_count = length(stat.time);
% number of temporal samples in the statistics object
j = [0:timestep:1]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];

%%
[i1,i2] = match_str(raweffectSLOWvsFAST.label, stat.label);

for k = 1:20;
   subplot(4,5,k);
   cfg = [];
   cfg.xlim=[j(k) j(k+1)];   % time interval of the subplot
   cfg.zlim = [-2.5e-13 2.5e-13];
   % If a channel reaches this significance, then
   % the element of pos_int with an index equal to that channel
   % number will be set to 1 (otherwise 0).

   % Next, check which channels are significant over the
   % entire time interval of interest.
   pos_int = zeros(numel(raweffectSLOWvsFAST.label),1);
   neg_int = zeros(numel(raweffectSLOWvsFAST.label),1);
   pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
   neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2);

   cfg.highlight = 'on';
   % Get the index of each significant channel
   cfg.highlightchannel = find(pos_int | neg_int);
   cfg.comment = 'xlim';
   cfg.commentpos = 'title';
   cfg.layout = '/mnt/home/a_shishkina/fieldtrip/template/layout/neuromag306planar.lay';
   cfg.interactive = 'no';
   ft_topoplotER(cfg, raweffectSLOWvsFAST);
end
