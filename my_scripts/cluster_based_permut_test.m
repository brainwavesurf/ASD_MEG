% ft_pms_step1_sensors.m
%
% Load data epochs 
% Load info about the preseeding trials
% Calculate spectral power on sensors in [-.8 to 0] interval
% plot average alpha power  distribution of epochs following different
% types of stimuli (Slow, Medium, Fast); plot Fast-Slow power differences
% save pictures as PPTX

%restoredefaultpath;
clear;
close all;
screensize = get( groot, 'Screensize' );

tapsmofrq  = 2;
megfolder = strcat( 'meg_sensors_tapsmofrq', num2str(tapsmofrq));
alpharange = [7 14];

%%
% NB: add to path: 
fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(path, fieldtripfolder)
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);

realdatapath = '/home/a_shishkina/data/KI/SUBJECTS/';
savepath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';
%%
%add list of subjects:
%SUBJ = ['0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0136'; '0137'; '0138'; '0139'; '0140'; '0141'; '0158'; '0159'; '0160'; '0161'; '0162'; '0163'; '0164'; '0178'; '0179'; '0253'; '0254'; '0255'; '0256'; '0257'; '0259'; '0273'; '0274'; '0275'; '0276'; '0277'; '0346'; '0347'; '0348'; '0350'; '0351'; '0357'; '0358'; '0378'; '0380'; '0381'; '0382'; '0383'; '0384'; '0385']; 
SUBJ = ['0076'];
s=1;
subj = SUBJ (s,:); 
savemegto = strcat(savepath, subj);
epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
%%
load ([ savemegto, '/', subj, '_info.mat'])
ev1 = find(allinfo.prev_stim_type==2); % i.e. epochs following Slow (2)
%ev2 = find(allinfo.prev_stim_type==4); %... medium
ev2 = find(allinfo.prev_stim_type==8); % ... fast
prev = {ev1, ev2};

ep_fiff_file = strcat(epofolder, subj, '-noerror-lagcorrected-epo.fif')
hdr = ft_read_header(ep_fiff_file);
%% Do the trial definition for slow visual grating condition:
fiff_file = strcat(realdatapath, subj, '/ICA_nonotch_crop/', subj, '_rings_raw.fif');
hdrraw = ft_read_header(fiff_file);
first = hdrraw.orig.raw.first_samp; %first sample length
hdr.orig.epochs.events(:,1) = hdr.orig.epochs.events(:,1) - first;
evraw = hdr.orig.epochs.events
%extract data about data from hdr
%evraw = [evraw(1,1); evraw([2:227],1) - evraw([1:226],1)]
%form trl with data corresponding to the index for correcevrawt previous stimulus
trlraw = evraw(prev{1},:)
%%
pre = -0.8*hdr.Fs;
post = -0.01*hdr.Fs;
trl = [];
for i=1:size (prev{1},2)
    trl(i, 1) = trlraw(i,1)+pre; 
    trl(i, 2) = trlraw(i,1)+post; 
    trl(i, 3) = 1.0*hdr.Fs ; % offset
     % stimulus_value;
end
%% Load unfiltered epochs. If you need to filter the data (e.g. for LCMV), import raw, not epochs.
ep_fiff_file = strcat(epofolder, subj, '-noerror-lagcorrected-epo.fif')
hdr = ft_read_header(ep_fiff_file);


cfg = [];  
cfg.dataset = ep_fiff_file;
cfg.channel={'megplanar'};
epochs = ft_preprocessing(cfg);

for con=1:2 % for conditions
        %select epochs according to preceding trials
        cfg = [];
        cfg.trials = prev{con}; % EV{con}';
        cfg.latency = [-0.8 0.0];
        dataPre{con} = ft_selectdata(cfg, epochs);
        
end    

dataSlow = dataPre{1};
dataFast = dataPre{2};
save dataSlow dataSlow
save dataFast dataFast

cfg = [];
cfg.planarmethod = 'sincos';
cfg.templates = 'neuromag306planar_neighb.mat';
% prepare_neighbours determines with what sensors the planar gradient is computed
cfg.method = 'distance';
neighbours = ft_prepare_neighbours(cfg, dataFast);
%%
dataSlow_planar = ft_megplanar(cfg, dataSlow);
dataFast_planar  = ft_megplanar(cfg, dataPre{2});
%%
cfg = [];
cfg.method = 'mtmfft';
cfg.output ='pow'; % 'fourier'
cfg.taper = 'hanning';
cfg.keeptrials = 'yes';
cfg.pad = 10; % to [vertually] increase spectral resolution
cfg.tapsmofrq = tapsmofrq;
cfg.foilim = alpharange; %freq band of interest
freqFast = ft_freqanalysis(cfg, dataFast); 
freqSlow = ft_freqanalysis(cfg, dataSlow);

%%
cfg = [];
freqFast_planar_cmb = ft_combineplanar(cfg, freqFast);
freqSlow_planar_cmb = ft_combineplanar(cfg, freqSlow);

%%
freqFast_planar_cmb.grad = dataFast.grad;
freqSlow_planar_cmb.grad = dataSlow.grad;

%%
cfg = [];
cfg.frequency = alpharange;
cfg.channel = {'MEG'}
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.computecritval = 'yes';
cfg.correctm = 'cluster';
cfg.clusteralpha  = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 500;
cfg.neighbours = neighbours;

design = zeros(1,size(freqSlow_planar_cmb.powspctrm,1) + size(freqFast_planar_cmb.powspctrm,1));
design(1,1:size(freqSlow_planar_cmb.powspctrm,1)) = 1;
design(1,(size(freqSlow_planar_cmb.powspctrm,1)+1):(size(freqSlow_planar_cmb.powspctrm,1)+...
size(freqFast_planar_cmb.powspctrm,1))) = 2;

cfg.design           = design;
cfg.ivar             = 1;

[stat] = ft_freqstatistics(cfg, freqSlow_planar_cmb, freqFast_planar_cmb);
%%
cfg = [];
freqFast = ft_freqdescriptives(cfg, freqFast_planar_cmb);
freqSlow = ft_freqdescriptives(cfg, freqSlow_planar_cmb);
%%
stat.raweffect = freqSlow.powspctrm - freqFast.powspctrm;

cfg = [];
cfg.alpha  = 0.025;
cfg.parameter = 'raweffect';
cfg.zlim   = [-1e-27 1e-27];
cfg.layout = 'neuromag306planar_helmet.mat';
ft_clusterplot(cfg, stat);


