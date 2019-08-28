% ft_pms_step1_sensors_and_DICS_permstat.m
% Do TMF on sensors and save sensor-level parameters in subject's folder
% (wF, p, etc.)
% Do source analysis for wF and save results in subject's folder
% Do clustervise stat to reveal significant gamma increases in the whole
% brain, as well as in occipital locations

% Uses precalculated  by  ft_step0_mri leadfields

%restoredefaultpath;
clear;
close all;
clc;

tapsmofrq  = 5;
gridres = 6; % 6 mm grid
megfolder = strcat( 'meg', num2str(gridres), 'mm_nonlinwarp_tapsmofrq', num2str(tapsmofrq), '_DICS');
mrifolder = 'mri_nonlinwarp_6mm_brthr0.5';
gammarange = [35 110];
lambda=5;
Nrand = 10000; % number of randomizations for bootstraping
Nvoxels = 25; % number of voxels [closest to the max] to average.

%%
% NB: add to path: 
fieldtripfolder = '/home/kolai/Documents/Shishkina/ProgramFiles/Fieldtrip/';
path(path, fieldtripfolder)
ft_defaults;
path('/home/kolai/Documents/Shishkina/ProgramFiles/Fieldtrip/external/mne/', path);


% This is an external matlab package used to save figures to PPTX
ft_defaults;
path('home/kolai/Documents/Shishkina/External_matlab/exportToPPTX_master/', path);

% EEEG lab functions
path (path, '/home/kolai/Documents/Shishkina/External_matlab/eeglab2019_0/functions/sigprocfunc');
path(path,'/home/kolai/Documents/Shishkina/External_matlab/eeglab2019_0/functions/popfunc');
path(path,'/home/kolai/Documents/Shishkina/External_matlab/eeglab2019_0/functions/guifunc');
path(path,'/home/kolai/Documents/Shishkina/External_matlab/eeglab2019_0/functions/adminfunc');
file = '/home/kolai/Documents/Shishkina/NeuralDataAnalysis/Autism/Code/scripts/for_meg_locs1.set'; % loc file
EEG = pop_loadset(file);
realdatapath = '/home/kolai/Documents/Shishkina/NeuralDataAnalysis/Autism/0101/MEGdata/';
DATAPATH  = '/home/kolai/Documents/Shishkina/NeuralDataAnalysis/Autism/0101/MEGdata/FT_beamf/'; % ???

%%%% Read tepllates
%%
% templatemodel
templatedir = strcat (fieldtripfolder , 'template/sourcemodel');
temp = strcat ('standard_sourcemodel3d', num2str(gridres), 'mm.mat')
template = load(fullfile(templatedir, temp)); 
%%%% template mri
templatefile = strcat (fieldtripfolder, '/external/spm8/templates/T1.nii'); 
template_mri = ft_read_mri(templatefile);
% atlas
atlas = ft_read_atlas( strcat (fieldtripfolder, '/template/atlas/aal/ROI_MNI_V4.nii') ); 
atlas = ft_convert_units(atlas,'cm');% assure that atlas and template_grid are expressed in the %same units

%%
                              
% SUBJ = [ 'G001'; 'G001';  'G002'; 'G003'; 'G003'; 'G005';'G006';'G007'; 'G007';  'G013'; 'G014'; 'G015'; 'G015'; 'G016'; 'G016'; 'G017'; 'G017'; 'G020'; 'G020'; 'G021'; 'G021'; 'G022'; 'G022'; 'G023'; 'G023'; 'G025'; 'G025'; 'G026'; 'G026'; 'G027'; 'G027'; 'G028'; 'G028'; 'G030'; 'G030'];
% DATE = [ 'f';    'l';      'l';       'f';    'l';    'l';   'f';   'f';   'l'; 'l'; 'l'; 'f'; 'l'; 'f'; 'l'; 'f'; 'l';  'f'; 'l'; 'f'; 'l'; 'f'; 'l'; 'f'; 'l';'f'; 'l'; 'f'; 'l'; 'f'; 'l'; 'f'; 'l'; 'f'; 'l'];
SUBJ = ['0101'];
DATE = ['oct-6'];
%%
s=1;
close all
subj = SUBJ (s,:); date = DATE (s,:);
DATAfolder = strcat(realdatapath);
MRIfolder = strcat(DATAPATH, subj, '_', date, '/', mrifolder);
savemegto = strcat(DATAPATH, subj, '_', date, '/', megfolder);
mkdir(savemegto)

%start PPTX report
exportToPPTX('new');
PPTXname  = strcat(savemegto, '/', subj, '_DISC_Tstat_report.pptx');

% template grid/leadfield
load(strcat(MRIfolder,'/', subj, '_template_grid.mat'));

%% Load subject's MEG and MRI:
% MNI warped leadfield grid
load(strcat( MRIfolder, '/', subj, '_grid_MNI_lf.mat'));
% volume/headmodel
load(strcat( MRIfolder, '/', subj, '_individ_hdm_vol.mat' )); % subject's template
% load subject mri: mri_orig_realigned
load (strcat (  MRIfolder, '/',  subj,  '_mri_orig_realigned.mat'))

%%     %% read epochs
% %     cfg = [];
% %     cfg.dataset = ep_fiff_file;
% %     data1 = ft_preprocessing(cfg);
%     %%
%     C = strsplit(hdr.orig.epochs.drop_log,'], ');
%     find1=strcmp(C, '[[');
%     find2=strcmp(C, '[');
%     find3=strcmp(C, '[]]');
%     ind = find(find1+find2+find3);

%% 
SUBJ2 = ['0102'];
subj2 = SUBJ2(s,:)

%
ep_fiff_file = strcat(DATAfolder, subj2, '-noerror-lagcorrected-epo.fif');
hdr = ft_read_header(ep_fiff_file);

ev_name=[realdatapath, 'epochs1/', subj2, date, '_events.mat']
load(ev_name)

%ev_name=['/home/kolai/Documents/Shishkina/NeuralDataAnalysis/Autism/0101/0101_clean_events.mat']
%load(ev_name);
% load epochs
cfg = [];
cfg.dataset = ep_fiff_file;
cfg.channel={'MEG*'};
epochs = ft_preprocessing(cfg);
%  select epochs according to events
    
ev1 = find(events(:,3)==2);
ev2 = find(events(:,3)==4);
ev3 = find(events(:,3)==8);
EV = {ev1, ev2, ev3};

 %% Load raw data
%     fiff_file = strcat(realdatapath, subj, '/ICA_notch_', date, '/', subj, '_static_raw.fif');
%     hdrraw = ft_read_header(fiff_file);
%     first= round(cast(hdrraw.orig.raw.first_samp, 'double'));
%     events(:,1) =  events(:,1)-first;
%     
%     % Define trials
%     % trl:   start, end and offset (interval before the event) 
%     pre = -1.0* hdr.Fs ;
%     post = 1.2* hdr.Fs ;
%     trl=[];
%     for i=1:size (events,1)
%         trl(i, 1)=(events(i,1)+pre) ; 
%         trl(i, 2)=(events(i,1)+post) ; 
%         trl(i, 3)= -1.0*hdr.Fs ; % offset
%         trl(i, 4) = events(i,3); % stimulus_value;
%     end
% 
%     % extract data and epochs from the raw
%     cfg = [];
%     cfg.trl=trl;
%     cfg.channel     = 'meg';
%     % cfg.demean      = 'yes';
%     cfg.dftfilter   = 'yes';
%     cfg.dftfreq     = [50 100];
%     cfg.dataset = fiff_file;
%     cfg    = ft_definetrial(cfg);
%     epochs = ft_preprocessing(cfg);
%     
%% Browse epochs
% %     cfg=[]
% %     cfg.viewmode = 'vertical';
% %     cfg.plotevents = 'yes';
% %     cfg.blocksize =5;
% %     ft_databrowser(cfg, epochs);

%% Plot average ERF
cfg = [];
cfg.channel = epochs.label;
avg = ft_timelockanalysis(cfg,epochs);
%%
hh=figure;
ttt = find(avg.time>-0.3 & avg.time<0.5);
plot (avg.time(ttt), avg.avg(:, ttt));
title ('Sensor averages')
exportToPPTX('addslide'); % slide with distributions
exportToPPTX('addpicture', hh, 'Position', [0.5,0.5,9,6]);
close (hh)

%% Calculating power in sensors, find max weighted freq, freq range and P of gamma enhancement in this freq range.
hh=figure;
LineColor{1} = 'k'; LineColor{2} = 'b'; LineColor{3} = 'g'; 
LineColorMax{1} = '--k'; LineColorMax{2} = '--b'; LineColorMax{3} = '--g'; 
ax=gca;

Ch = {'MEG1932',  'MEG1922', 'MEG2042',  'MEG2032',  'MEG2112', 'MEG2122',  'MEG2342', 'MEG2332',  'MEG1732', 'MEG1942', 'MEG1912', 'MEG2012', 'MEG2022', 'MEG2312', 'MEG2322', 'MEG2512',...
      'MEG1933',  'MEG1923', 'MEG2043',  'MEG2033',  'MEG2113', 'MEG2123',  'MEG2343', 'MEG2333',  'MEG1733', 'MEG1943', 'MEG1913', 'MEG2013', 'MEG2023', 'MEG2313', 'MEG2323', 'MEG2513'};

exportToPPTX('addslide'); % slide with distributions
exportToPPTX('addtext','Power ratio [post/pre] distribution for GRA1 and GRA2, conditions 1,2,3');
exportToPPTX('addtext','GRA1 sensor set, ... MEG0112,...', 'Position', [0.5, 0.5, 10, 1]);
exportToPPTX('addtext','GRA2 sensor set, ... MEG0113,...', 'Position', [0.5, 3.6, 10, 1]);

X=[0.4, 2.6, 4.8, 7.0; 0.4, 2.6, 4.8, 7.0]; Y=[1, 1, 1, 1; 4, 4, 4, 4]; % positions for plot
%%
for con=1:3 % for conditions
    cfg = [];
    cfg.trials = EV{con}; % EV{con}';
    cfg.latency = [-0.8 0.0]
    [dataPre] = ft_selectdata(cfg, epochs);
    cfg.latency = [0.4 1.2]
    [dataPost] = ft_selectdata(cfg, epochs);

    cfg = [];
    cfg.method    = 'mtmfft';
    cfg.output    ='pow'; % 'fourier'
    cfg.taper        = 'dpss';
    cfg.keeptrials = 'yes';
    cfg.tapsmofrq = tapsmofrq;
    cfg.foilim    = gammarange; %freq band of interest
    freqPre = ft_freqanalysis(cfg, dataPre);  % trials x Ch x freq
    freqPost = ft_freqanalysis(cfg, dataPost);
    % N of max channel
    for j=1:length(Ch)
        [ch(j),x] = find(strcmp(freqPost.label,Ch{j}));
    end
    %% Find peak freq based on max probability value
    % b1 is the max ch
    [a,b1] = max(squeeze(mean(mean(freqPost.powspctrm(:,ch,:),3),1))./squeeze(mean(mean(freqPre.powspctrm(:,ch,:),3),1))); % trial x ch x freq
    ratio{con} =  (squeeze(mean(freqPost.powspctrm(:,ch(b1),:),1))./squeeze(mean(freqPre.powspctrm(:,ch(b1),:),1)) -1)'; % post/pre ratio at the 'max channel'
    ChName{con}  = freqPost.label(ch(b1));

    % find N=4 max channels
    change =squeeze(mean(freqPost.powspctrm(:,ch,:),1))./squeeze(mean(freqPre.powspctrm(:,ch,:),1)); % ch x freq;
    change = mean(change,2); % over frequencies in gamma range
    change =[1:length(change)'; change']';
    change = sortrows(change,2);
    selection4=change([end, end-1, end-2, end-3],1); %N channels with max (max-1, max-2, etc) power CHANGE in the GAMMA BAND
    t= ch (selection4);

    Ratio{con} =  (  squeeze(mean(mean(freqPost.powspctrm(:,t,:),1),2))./squeeze(mean(mean(freqPre.powspctrm(:,t,:),1),2) ) -1   )'; % post/pre ratio at the 'max channel selection'
    plot (ax,   freqPost.freq, Ratio{con}, LineColor{con}); hold on, plot(ax, freqPost.freq, ratio{con},  LineColorMax{con});

    bins{con} = [];
    for f=1:size(freqPost.powspctrm,3) % for all freq
         post = squeeze(freqPost.powspctrm(:,ch(b1),f)); % ch(b1) is the max ch
         pre  = squeeze(freqPre.powspctrm(:,ch(b1),f));    
         [p_f{con}(f),h,stats] = ranksum(pre,post);
    end
    [val, ind ] = min(p_f{con})
    PeakF_Pbased{con}(1) = freqPre.freq(ind);
    PeakF_Pbased{con}(2) = val;

    %% Find peak freq based on absolute max value

    %find max weighted freq
    [MAXratio, MAXind] = (max(Ratio{con})); maxF = freqPost.freq(MAXind);
    bins{con} = find( (Ratio{con}) > max(Ratio{con} )/3*2); % 'bins' is the subject/condition-specific band
    binsfreqs{con} = freqPost.freq(bins{con});

    % if frequencies are too far from the maximum (more than 20 Hz) we exclude them!
    excl = find(abs(binsfreqs{con}-maxF)>20);
    bins{con}(excl) =[];
    binsfreqs{con} = freqPost.freq(bins{con});

    wF{con} = sum(binsfreqs{con}.*Ratio{con}(bins{con}))/sum(Ratio{con}(bins{con}));  % if 'Ratio', then at the 'max channel selection, N=4, see line 171'
    wPmax{con} = mean(ratio{con}(bins{con})+1); % max ch
    wPsel{con} = mean(Ratio{con}(bins{con})+1); % selection of 4 channels

    % in the max channels find P of gamma increase
    [aa,bb] = min(abs(freqPost.freq-wF{con})); % bb is the closest max freq    

    % for subject/condition-specific band
    post = squeeze(mean(mean(freqPost.powspctrm(:,t,bins{con}),3),2));  %bins{con} t is the max ch selection (of 4 ch), at the subject/condition-specific band
    pre  = squeeze(mean(mean(freqPre.powspctrm(:,t,bins{con}),3) ,2));
    [p_band{con},h,stats] = ranksum(pre,post);
    postave{con} = mean(post)*1e24;
    preave{con}  = mean(pre)*1e24;

    % for maximum
    post = squeeze(mean(freqPost.powspctrm(:,t,bb),2));  %bb is the nearest to wF frequency bin
    pre = squeeze(mean(freqPre.powspctrm(:,t,bb),2));       
    [p_max{con},h,stats] = ranksum(pre,post);

%       postave{con} = squeeze(mean(mean(freqPost.powspctrm(:,(ch(b)),bins{con}),3),1))*1e24;
%       preave{con}  = squeeze(mean(mean(freqPre.powspctrm(:,(ch(b)),bins{con}),3),1))*1e24;
%%
    data = squeeze(mean(mean(freqPost.powspctrm(:,:,ind),3),1))./squeeze(mean(mean(freqPre.powspctrm(:,:,ind),3),1))-1; % ind is the p-based peak freq
    [gra,y] = find(strcmp(freqPost.grad.chantype, 'megplanar'));

    data = data(gra);
    Data=data(1:2:length(data));
    h=figure;
    topoplot (Data, EEG.chanlocs, 'electrodes', 'on', 'maplimits', [0,max(ratio{1})]);    colorbar
    title (strcat(subj, ',..GRA1..V', num2str(con), ',..maxCh=', ChName{con}),  'FontSize', 20)
    exportToPPTX('addpicture', h, 'Position', [X(1,con),Y(1,con),2.3,1.8]);
    Data=data(2:2:length(data));
    close (h)
    h=figure;
    topoplot (Data, EEG.chanlocs, 'electrodes', 'on','maplimits', [0,max(ratio{1})]);    colorbar
    title (strcat(subj, ',..GRA2..V', num2str(con), ',..maxCh=', ChName{con}), 'FontSize', 20)
    exportToPPTX('addpicture', h, 'Position', [X(2,con),Y(2,con),2.3,1.8]);
    close (h)
end  % end for conditions in sensor space
%%
title(ax, strcat('[post-pre]/pre ratio in MAX-prob and 3-AVE channels:',  ChName{1}, '../',ChName{2}, '../', ChName{3}),  'FontSize', 10 )
legend(ax, {'SlowAVE', 'SlowMAX', 'Medium','MediumMAX',  'FastAVE', 'FastMAX'}, 'FontSize', 10);
exportToPPTX('addslide'); % slide with distributions
exportToPPTX('addpicture', hh, 'Position', [1,1,6,4]);
exportToPPTX('addtext', strcat('F1=', num2str(wF{1}), ', F2=', num2str(wF{2}), ', F3=', num2str(wF{3})), 'FontSize', 11, 'Position', [0.5, 5.0, 10, 1]);
exportToPPTX('addtext', strcat('P1max_dB=', num2str(wPmax{1}),'P2max_dB=', num2str(wPmax{2}), ', P3max_dB=', num2str(wPmax{3})), 'FontSize', 11, 'Position', [0.5, 5.4, 10, 1]);
exportToPPTX('addtext', strcat('P1selection4_dB=', num2str(wPsel{1}), ', P2selection4_dB=', num2str(wPsel{2}), ', P3selection4_dB=', num2str(wPsel{3})), 'FontSize', 11, 'Position', [0.5, 5.6, 10, 1]);

exportToPPTX('addtext', strcat('Condition-specific band for V1:...', num2str(binsfreqs{1}(1)), '..to..', num2str(binsfreqs{1}(length(binsfreqs{1})))   ), 'FontSize', 11,  'Position', [0.5, 6.0, 10, 1]);
exportToPPTX('addtext', strcat('Condition-specific band for V2:...', num2str(binsfreqs{2}(1)), '..to..', num2str(binsfreqs{2}(length(binsfreqs{2})))   ), 'FontSize', 11,  'Position', [0.5, 6.2, 10, 1]);
exportToPPTX('addtext', strcat('Condition-specific band for V3:...', num2str(binsfreqs{3}(1)), '..to..', num2str(binsfreqs{3}(length(binsfreqs{3})))   ), 'FontSize', 11,  'Position', [0.5, 6.4, 10, 1]);

exportToPPTX('addtext', strcat('M-W rank test at freq closest to wF:  pV1_max=', num2str(p_max{1}), ', pV2_max=', num2str(p_max{2}), ', pV3_max=', num2str(p_max{3})), 'FontSize', 11,  'Position', [0.5, 7.1, 10,1]);
exportToPPTX('addtext', strcat('M-W rank test in subj/freeq-specific band:  pV1_band=', num2str(p_band{1}), ', pV2_band=', num2str(p_band{2}), ', pV3_band=', num2str(p_band{3})), 'FontSize', 11, 'Position', [0.5, 7.3, 10,1]);
close (gcf)


%% Save sensor info
readme = [];
readme.wF = 'Individual weighterd frequency, at 2/3 max in max ch';
readme.wPmax = 'Max weighted pow in dB in max channel ';
readme.wPsel = 'Max weighted pow in dB in the selection-of-4-max-channels';
readme.p_band = 'probability at the selection-of-4-max-channels amplitude increase ';
readme.p_max = 'probability at the selection-of-4-max-channels amplitude increase, is the nearest to wF frequency bin ';
readme.ChName = 'Max sensor name';
readme.preave = 'prestim power in max sensor selection, averaged over condition-specific band';
readme.postave = 'poststim power in max sensor selection, averaged over condition-specific band';
readme.binsfreqs = 'condition specific bands';
readme.PeakF_Pbased= 'peak F at max probability and this probability';

Savename = strcat(savemegto,'/', subj, '_sensors_res.mat');
% Max weighted freq
save (Savename, 'wF');
% Max weighted pow in dB in max channel
save (Savename, 'wPmax', '-append');
% Max weighted pow in dB in the channel selection (n=4)
save (Savename, 'wPsel', '-append');
% P value for condition-specific band
save (Savename, 'p_band', '-append');
% P value for the nearest max freq
save (Savename, 'p_max', '-append');
% Max Sensor
save (Savename, 'ChName', '-append');
% prestim power in max sensor selection, averaged over condition-specific band
save (Savename, 'preave', '-append'); 
% poststim power in max sensor selection, averaged over condition-specific band
save (Savename, 'postave', '-append'); 
% condition specific bands
save (Savename, 'binsfreqs', '-append');   
save (Savename, 'PeakF_Pbased', '-append');   

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   sourse localization of MEG data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating the cross spectral density matrix (The DICS spatial filter is derived from the frequency counterpart of the covariance matrix)

for con= 1:3  % for conditions
    cfg = [];
    cfg.trials = EV{con}; % EV{con}';
    cfg.latency = [-0.8 0.0]
    [dataPre] = ft_selectdata(cfg, epochs);
    cfg.latency = [0.4 1.2]
    [dataPost] = ft_selectdata(cfg, epochs);

    %% critical t value
    ntrials = numel(dataPost.trial);
    Tval = tinv(1-0.025,ntrials-1);

    %% spectral analysis, sensors
    cfg = [];
    cfg.method    = 'mtmfft';
    cfg.output    = 'fourier'; %'powandcsd'; %
    cfg.keeptrials = 'yes';
    cfg.grad         = dataPre.grad;
    cfg.tapsmofrq = tapsmofrq;
    cfg.foi     = wF{con}; %freq band of interest
    %cfg.foi     = 70;
    freqPre = ft_freqanalysis(cfg, dataPre);
    freqPost = ft_freqanalysis(cfg, dataPost);
    dataAll = ft_appenddata([], dataPre, dataPost);
    freqAll = ft_freqanalysis(cfg, dataAll);

    %% Then we compute the inverse filter based on both [pre vs post] conditions. 
    %  Note the use of cfg.keepfilter so that the output saves this computed filter.
    cfg              = [];
    cfg.method       = 'dics';
    cfg.frequency    = wF{con};
    cfg.grid         = grid_MNI_lf; %fwd model,  see FT_PREPARE_SOURCEMODEL or FT_PREPARE_LEADFIELD
    cfg.headmodel    = individ_hdm_vol;
    cfg.dics.lambda       = lambda;
    cfg.dics.fixedori  = 'yes'; %is that we only keep the largest of the three dipole directions per spatial filter
    cfg.dics.keepfilter   = 'yes'; % we let ft_sourceanalysis return the filter matrix in the source structure. 
    cfg.dics.realfilter   = 'yes'; % specifies that we do not allow our filter to have an imaginary part. 
    sourceAll = ft_sourceanalysis(cfg, freqAll);
    %sourceAll.pos = template.sourcemodel.pos;  %template.pos;

% %         %% tmp
% %         %  To plot an 'orthogonal cut' in subject's coordinates:
% %         cfg = [];
% %         cfg.method        = 'ortho';
% %         cfg.funparameter  = 'avg.pow';
% %         cfg.maskparameter = cfg.funparameter;
% %         cfg.funcolorlim   = [-1*max(sourceAll.avg.pow) max(sourceAll.avg.pow)];
% %         %cfg.opacitylim    = [Tval Tval*2]; 
% %         cfg.opacitymap    = 'rampup'; 
% %         ft_sourceplot(cfg, sourceAll);   

    %%  By placing this pre-computed filter inside cfg.grid.filter, it can now be applied to each condition separately.             
    cfg=[];
    cfg.method       = 'dics';
    cfg.frequency    = wF{con};
    cfg.grid         =  grid_MNI_lf; %fwd model,  see FT_PREPARE_SOURCEMODEL or FT_PREPARE_LEADFIELD
    cfg.headmodel    = individ_hdm_vol;
    cfg.grid.filter = sourceAll.avg.filter;
    cfg.fixedori  = 'yes';
    cfg.rawtrial = 'yes';
    % Here we do it for single trials to perform statistics later on:
    sourcePre_con  = ft_sourceanalysis(cfg, freqPre );  %sourcePre_con.avg.pow is pow: [11000, 1 double]
    sourcePost_con = ft_sourceanalysis(cfg, freqPost);

    % Here we do the same, but for average power, and later we save the ave power:
    cfg.rawtrial = 'no';
    sourcePre_ave_con  = ft_sourceanalysis(cfg, freqPre );  %sourcePre_con.avg.pow is pow: [11000?1 double], i.e. average in the band ?!
    sourcePost_ave_con = ft_sourceanalysis(cfg, freqPost);
    filename = strcat(savemegto, '/', subj, strcat('_sourcePre_wF_ave_con', num2str(con), '.mat')  );
    save (filename, 'sourcePre_ave_con')
    filename = strcat(savemegto, '/', subj, strcat('_sourcePost_wF_ave_con', num2str(con), '.mat')  );
    save (filename, 'sourcePost_ave_con')

    %%  Now we can compute the contrast of (post-pre)/pre. 
% %         % In this operation we assume that the noise bias is the same for the pre- and post-stimulus interval and it will thus be removed.
% %         sourceDiff = sourcePost_con; % [3276x3 double]
% %         sourceDiff.avg.pow = (sourcePost_con.avg.pow - sourcePre_con.avg.pow) ./ sourcePre_con.avg.pow;

%% Now we can statistically compare the difference between the pre and post response source power in the WHOLE BRAIN using ft_sourcestatistics 
cfg = [];
cfg.parameter    = 'pow';
cfg.dim          = sourcePost_con.dim;
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';  %'ft_statfun_depsamplesT'; %
cfg.correctm         = 'max'; %'cluster'; %
%      cfg.clusterstatistic = 'maxsum';
%      cfg.clusteralpha     = 0.05;
%      cfg.clusterthreshold = 'nonparametric_individual';
%      cfg.clustercritval = 1;% for parametric thresholding (default is determined by the statfun)
%      cfg.clustertail      = 1; % (i.e. our hypothesis is that post > pre)
cfg.tail             = 1; % (i.e. our hypothesis is that post > pre)
cfg.alpha            = 0.05; % (i.e. our hypothesis is that post > pre)
cfg.numrandomization = Nrand;

design  = zeros(2,2*ntrials);
design(1,1:ntrials) = 1;
design(1,ntrials+1:2*ntrials) = 2;
design(2,1:ntrials) = 1:ntrials;
design(2,ntrials+1:2*ntrials) = 1:ntrials;

cfg.design   = design;
cfg.ivar     = 1;
cfg.uvar     = 2;
stat = ft_sourcestatistics(cfg,sourcePost_con,sourcePre_con);
stat_templ = stat;
stat.pos=template_grid.pos;% keep positions for plotting later
%stat.unit='cm';
STAT{con}= stat;

% %     %% tmp
% %     %  To plot an 'orthogonal cut' in subject's coordinates:
% %     cfg = [];
% %     cfg.method        = 'ortho';
% %     cfg.funparameter  = 'stat';
% %     cfg.maskparameter = cfg.funparameter;
% %     cfg.funcolorlim   = [-Tval*2 Tval*2];
% %     %cfg.opacitylim    = [Tval Tval*2]; 
% %     cfg.opacitymap    = 'rampup'; 
% %     ft_sourceplot(cfg, stat_templ);

%% Interpolate atlas to functional (we do it one time)
if con==1
   cfg              = [];
   cfg.voxelcoord   = 'no';
   cfg.parameter    = 'tissue';
   cfg.interpmethod = 'nearest';
   AtlasInt  = ft_sourceinterpolate(cfg, atlas, stat);   
   AtlasInt.coordsys = 'mni';

   % By doing this you get the label for every grid  :
   % AtlasInt is the atlas interpolated to our functional data

   AtlasInt.tissue(isnan(AtlasInt.tissue)) =0; % replace NAN with 0 ;
   ids  =  find(AtlasInt.tissue);          % all interpolate regions
   id   =  AtlasInt.tissue(ids); %  all interpolate regions indexes
   ROI     = AtlasInt.tissuelabel(id);
    occid1   = find(strncmpi(ROI,'Occipital',9));  %  indice
    occid2   = find(strncmpi(ROI,'Calcarine',9));  %  indice
    occid3   = find(strncmpi(ROI,'Cuneus',6));  %  indice
    occid4   = find(strncmpi(ROI,'Lingual',7));  %  indice
    occid5   = find(strncmpi(ROI,'Precuneus', 9));  
    occid    = sort([occid1, occid2, occid3, occid4, occid5]);        
    %OCC      = ROI(occid);  % label
end

%% And now the same in the OCCIPITAL  mask
sourcePost_con_occ = sourcePost_con ;
sourcePost_con_occ.inside = zeros(size(sourcePost_con.inside));
sourcePost_con_occ.inside(ids(occid))  = 1;

sourcePre_con_occ = sourcePre_con ;
sourcePre_con_occ.inside = zeros(size(sourcePre_con.inside));
sourcePre_con_occ.inside(ids(occid))  =1; 

cfg = [];
cfg.parameter    = 'pow';
cfg.dim          = sourcePost_con_occ.dim;
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';  %'ft_statfun_depsamplesT'; %
cfg.correctm         = 'max'; %'cluster';
%     cfg.clusterstatistic = 'maxsum';
%     cfg.clusteralpha     = 0.05;
%     cfg.clustercritval = 1;% for parametric thresholding (default is determined by the statfun)
%     cfg.clustertail      = 1;
cfg.tail             = 1;
cfg.alpha            = 0.05;
cfg.numrandomization = Nrand;

cfg.design   = design;
cfg.ivar     = 1;
cfg.uvar     = 2;
stat_occ = ft_sourcestatistics(cfg,sourcePost_con_occ,sourcePre_con_occ);
stat_occ.pos=template_grid.pos;%template.sourcemodel.pos; %template_grid.pos;% keep positions for plotting later
STAT_OCC{con}= stat_occ;

%% Make the low resolution MRI corresponding to the functional data 
anatomy = ft_convert_units(template_mri, 'cm');
cfg              = [];
cfg.parameter    = 'anatomy';
%cfg.voxelcoord   = 'no';
cfg.interpmethod = 'nearest';
anatomy = ft_sourceinterpolate(cfg,  anatomy, stat);  
anatomy.coordsys = 'mni';

%     %%  tmp
%         %  To plot an 'orthogonal cut' in subject's coordinates:
%         cfg = [];
%         cfg.method        = 'ortho';
%         cfg.funparameter  = 'stat';
%         cfg.maskparameter = cfg.funparameter;
%         cfg.funcolorlim   = [-Tval Tval*2];
%         cfg.opacitylim    = [-Tval Tval*2]; 
%         cfg.opacitymap    = 'rampup'; 
%         ft_sourceplot(cfg, stat);   
%      
%% Subsequently we interpolate the result and the binary mask containing information of significant deferences per voxel.
 if size(find(stat_templ.mask==1),1)>0
    Cluster(con)=1;
    cfg              = [];
    cfg.voxelcoord   = 'no';
    cfg.parameter    = 'stat';
    cfg.interpmethod = 'nearest';
    statint  = ft_sourceinterpolate(cfg, stat, anatomy); %interpolate to anatomy
    cfg.parameter    = 'mask';
    maskint  = ft_sourceinterpolate(cfg, stat, anatomy);
    statint.mask = maskint.mask;

    % And plot the result masked for significant activations, the functional data is now expressed in t-values.
    statint.coordsys = 'mni';
    cfg               = [];
    cfg.method        = 'ortho';
    cfg.funparameter  = 'stat';
    cfg.maskparameter = 'mask';
    cfg.atlas         =  AtlasInt; % atlas;%
    cfg.location = 'max';
    cfg.funcolorlim   = [Tval Tval*2];
    cfg.funcolormap = 'hot';
    ft_sourceplot(cfg,statint);

    title (strcat(subj, '...V', num2str(con), '...cluster in the whole brain'), 'FontSize', 12);
    exportToPPTX('addslide');
    exportToPPTX('addtext', strcat( 'Condition...V', num2str(con), ', Source localization of T-stat in the whole brain, at central freq=', num2str(wF{con}), 'Hz '));
    exportToPPTX('addpicture', gcf,  'Position', [1,1,7.5,5]);
 else
    Cluster(con)=0;
    exportToPPTX('addslide');
    exportToPPTX('addtext', strcat( 'Condition...V', num2str(con), ', Whole brain correction, method=max: no significant voxels, at central freq=', num2str(wF{con}), 'Hz '));
 end
 close (gcf)

%% The same for occipital: Subsequently we interpolate the result and the binary mask containing information of significant deferences per voxel.
 if size(find(stat_occ.mask==1),1)>0
    Occ_Cluster(con)=1;
    cfg              = [];
    cfg.voxelcoord   = 'no';
    cfg.parameter    = 'stat';
    cfg.interpmethod = 'nearest';
    statint  = ft_sourceinterpolate(cfg, stat_occ, anatomy); %interpolate to 
    cfg.parameter    = 'mask';
    maskint  = ft_sourceinterpolate(cfg, stat_occ, anatomy);
    statint.mask = maskint.mask;

    % And plot the result masked for significant activations, the functional data is now expressed in t-values.
    statint.coordsys = 'mni';
    cfg               = [];
    cfg.method        = 'ortho';
    cfg.funparameter  = 'stat';
    cfg.maskparameter = 'mask';
    cfg.atlas         =  AtlasInt; % atlas;%
    cfg.location = 'max';
    cfg.funcolorlim   = [Tval Tval*2];
%         cfg.funcolormap = 'hot';
    ft_sourceplot(cfg,statint);

    title (strcat(subj, '...V', num2str(con), '...cluster in occipital'), 'FontSize', 12);
    exportToPPTX('addslide');
    exportToPPTX('addtext', strcat( 'Condition...V', num2str(con), ', Source localization of T-stat in the occipital region, at freq of weighted max at...', num2str(wF{con}), 'Hz '));
    exportToPPTX('addpicture', gcf,  'Position', [1,1,7.5,5]);
 else
    Occ_Cluster(con)=0;
    exportToPPTX('addslide');
    exportToPPTX('addtext', strcat( 'Condition...V', num2str(con), ',  Occipital ROIs correction, method=max: no significant voxels, at central freq=', num2str(wF{con}), 'Hz '));
 end
 close (gcf)

%% When combining the source-level estimates of activity in multiple subjects, the activity can first be interpolated on the individuals MRI 
% (using ft_sourceinterpolate) 
cfg              = [];
%     cfg.voxelcoord   = 'no';
cfg.parameter    = 'stat';
cfg.interpmethod = 'nearest';
sourceStatInt  = ft_sourceinterpolate(cfg, stat_templ, mri_orig_realigned );   %mri_orig mri_orig_realigned

% ... and then spatially normalized to a template brain (using ft_volumenormalise). 
tmp =  ft_convert_units(sourceStatInt, 'mm');
cfg = [];
cfg.nonlinear     = 'yes';
sourceStatIntNorm = ft_volumenormalise(cfg, tmp);
sourceStatIntNorm = ft_convert_units(sourceStatIntNorm, 'cm');

% Here we interpolate to template
cfg              = [];
cfg.parameter    = 'stat';
cfg.interpmethod = 'nearest';
sourceStatIntTempl  = ft_sourceinterpolate(cfg, stat, template_mri);   
sourceStatIntTempl = ft_convert_units(sourceStatIntNorm, 'cm');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%
%  Plotting

%% Now plot T values using 'slice method':
% Slices in subject[s coordinates:
cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'stat';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [-Tval*2 Tval*2];
%cfg.opacitylim    = [Tval Tval*2]; 
cfg.opacitymap    = 'rampup';  
ft_sourceplot(cfg, sourceStatInt);
title (strcat(subj, '...V', num2str(con), '...sourceplot'), 'FontSize', 17);
exportToPPTX('addslide');
exportToPPTX('addtext', strcat( 'Condition...V', num2str(con), ', Source localization in weighted max at...', num2str(wF{con}), 'Hz '));
exportToPPTX('addpicture', gcf,  'Position', [-0.25,1,4,3]);
close (gcf)

%%  To plot an 'orthogonal cut' in subject's coordinates:
cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'stat';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [-Tval*2 Tval*2];
%cfg.opacitylim    = [Tval Tval*2]; 
cfg.opacitymap    = 'rampup'; 
ft_sourceplot(cfg, sourceStatInt);
exportToPPTX('addpicture', gcf,  'Position', [-0.1, 4.5, 3.5, 2.8]);
close (gcf)

%% The same in SPM coordinates
cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'stat';
cfg.funcolorlim   = [Tval Tval*2];
cfg.opacitylim    = [Tval Tval*2];
cfg.opacitymap    = 'rampup'; 
cfg.atlas         = atlas;
ft_sourceplot(cfg, sourceStatIntTempl);
title ('Normalized, Orthogonal view with  atlas')
exportToPPTX('addpicture', gcf,  'Position', [3.2 ,4.5, 3.5, 2.8]);
close (gcf)

%% .. and then on surface
    % sourceplot with method ?ortho? aftccccer volume normalisation
    % You can also project the power onto a surface using ft_sourceplot. FieldTrip has several surface .mat files available. 
    %The surface files are in MNI coordinates, so therefore the volume has to be normalized to match those coordinates. 
    %This can be done with the FieldTrip function ft_volumenormalise (see above, as well).

%        % NB : but we already normalized, because we used template grid!!!
    %normalises anatomical and functional volume data  to a template anatomical MRI.
    %[sourceDiffIntNorm] = ft_convert_units(sourceDiffIntNorm,'mm');
    sourceStatIntNorm1 = rmfield(sourceStatIntNorm,'coordsys'); % to some strange reason it does not work othervise!
    cfg = [];
    cfg.method         = 'surface';
    cfg.funparameter   = 'stat';
    cfg.maskparameter  = cfg.funparameter;
    cfg.funcolorlim    = [Tval Tval*2];
    cfg.funcolormap    = 'hot';
    cfg.opacitylim     = [Tval Tval*2]; 
    cfg.opacitymap     = 'rampup';  
    cfg.projmethod     =  'nearest'; %'sphere_avg'; 'project';
    cfg.surffile       = strcat(fieldtripfolder, '/template/anatomy/surface_white_both.mat'); % surface_inflated_both.mat surface_white_both triangulation that corresponds with SPM anatomical template in MNI coordinates
    cfg.surfdownsample = 10; 
    cfg.camlight ='no';
    %cfg.projthresh = 0.5
    cfg.sphereradius = 10
    ft_sourceplot(cfg, sourceStatIntNorm1); %
    view ([34 -14])

    title (strcat(subj, '...condition','...1', '...surface..sourceplot'));
    exportToPPTX('addpicture', gcf,  'Position', [3.33,0.5,4,3]);
    close (gcf)

    ft_sourceplot(cfg,  sourceStatIntNorm1); %
    view ([-34 -14])
    title (strcat(subj, '...condition','...1', '...surface..sourceplot'));
    exportToPPTX('addpicture', gcf,  'Position', [6.16,0.5,4,3]);
    close (gcf)

%% Interpolate to template mri
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'stat';
cfg.interpmethod = 'nearest';
sourceStatIntTemp  = ft_sourceinterpolate(cfg, stat, template_mri);
sourceStatIntTemp.coordsys= 'mni';
sourceStatIntTemp=  ft_convert_units(sourceStatIntTemp,'cm')

cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'stat';
cfg.funcolormap = 'hot';
cfg.funcolorlim   = [Tval Tval*2];
cfg.opacitylim    = [Tval Tval*2];
cfg.opacitymap    = 'rampup'; 
cfg.atlas         = atlas;
ft_sourceplot(cfg, sourceStatIntTemp);
title ('Interpolated to template, Orthogonal view with  atlas')
exportToPPTX('addpicture', gcf,  'Position', [3.2 ,4.5, 3.5, 2.8]);
close (gcf)

 %% Parcellate% 
cfg=[];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'stat';
cfg.interpmethod = 'nearest';
sourceStatIntTemp  = ft_sourceinterpolate(cfg, sourceStatIntTemp,  atlas);
parcel = ft_sourceparcellate(cfg, sourceStatIntTemp, atlas);
PARCEL{con}=parcel;  % save parcellation averages

%% Find max T in the whole brain and in OCC, record this information to be saved later on
[ maxT, maxind]= max(stat.stat(ids))
MaxMNICoord{con} = stat.pos(ids(maxind),:);
MaxMNI_Tstat{con} = maxT; 
if AtlasInt.tissue(ids(maxind))==0
   MaxisIn{con} = 'No_Label';
else
   a=AtlasInt.tissuelabel(AtlasInt.tissue(ids(maxind)));
   MaxisIn{con} = a{1};
end
POW_MaxMNI_Tstat{con}  = sourcePre_ave_con.avg.pow(ids(maxind));
POWpre_MaxMNI{con}  = sourcePre_ave_con.avg.pow(ids(maxind));
POWpost_MaxMNI{con}  = sourcePost_ave_con.avg.pow(ids(maxind));

[maxstatOCC, tmpind]= max(stat.stat(ids(occid)));
maxindOCC{con}(1) = ids(occid(tmpind));  maxindOCC{con}(2) = occid(tmpind);

MaxMNICoordOCC{con} = stat.pos(ids(occid(tmpind)),:);
MaxMNI_TstatOCC{con} = maxstatOCC; 
POWpre_MaxMNI_OCC{con}  = sourcePre_ave_con.avg.pow(ids(occid(tmpind)));
POWpost_MaxMNI_OCC{con}  = sourcePost_ave_con.avg.pow(ids(occid(tmpind)));

a= AtlasInt.tissuelabel(AtlasInt.tissue(ids(occid(tmpind))));
MaxisInOCC{con} =a{1};


%% Nvoxels average
ins = find(template_grid.inside);
x1=template_grid.pos(ins,1); y1=template_grid.pos(ins,2);  z1=template_grid.pos(ins,3); 
x0=template_grid.pos(occid(tmpind),1); y0=template_grid.pos(occid(tmpind),2);  z0=template_grid.pos(occid(tmpind),3); 
distance = sqrt(  (x1-x0).^2 +(y1-y0).^2 +(z1-z0).^2 ); 
distance = [ins , distance];
Cocc=sortrows(distance,2);
Nvoxels_id =Cocc(1:Nvoxels,:); Nvoxels_id = Nvoxels_id(:,1); Nvoxels_id = ins(Nvoxels_id );

POWpre_selMNI_OCC{con}  = sourcePre_ave_con.avg.pow(Nvoxels_id);
POWpost_selMNI_OCC{con}  = sourcePost_ave_con.avg.pow(Nvoxels_id);


    %% and plot parcellation
% We create a dummy structure where we identify the stat T values per voxel and use this for subsequent plotting. 
source_norm = sourceStatIntTemp ;
source_norm.anatomy = template_mri.anatomy
dummy=atlas;
for i=1:length(parcel.stat)
      dummy.tissue(find(dummy.tissue==i))=parcel.stat(i);
end;
source_norm.parcel=dummy.tissue;
source_norm.coordsys = 'mni';
cfg=[]; 
cfg.method = 'ortho';
cfg.funparameter = 'parcel';
cfg.funcolormap    = 'jet';
cfg.renderer = 'zbuffer';
cfg.location = [-0.2 -10 0];
cfg.atlas = atlas;
%cfg.funcolorlim = [-1 1];
ft_sourceplot(cfg,source_norm );  
title (strcat(subj, ', parcellated'));
exportToPPTX('addpicture', gcf,  'Position', [6.3, 4.5, 3.5, 2.8]);
close (gcf)

end  % end for conditions in source space

%% now save Source info:
filename = strcat(savemegto, '/', subj, '_DISCbroad_source_maxT.mat');
readme = [];
readme.MaxMNICoord = 'Coordinate of the maximal voxel in MNI space';
readme.MaxMNI_Tstat = 'Max T value';
readme.MaxisIn = 'Label of the max voxel';
readme.MaxMNICoordOCC = 'Coordinate of the maximal voxel in the occpital/visual regions, in MNI space';
readme.MaxMNI_TstatOCC = 'max T value, in the occpital/visual regions';
readme.MaxisInOCC = 'Label of the max voxel, in the occpital/visual regions';
readme.maxindOCC = 'index of the max voxel [1.in the whole template; 2. inside], in the occpital/visual regions'
readme.POWpre_MaxMNI_OCC = 'power at the max T value in occipital, in prestimulus'
readme.POWpost_MaxMNI_OCC = 'power at the max T value in occipital, in poststimulus'
readme.POWpre_MaxMNI = 'power at the max T value, in prestimulus'
readme.POWpost_MaxMNI = 'power at the max T value, in poststimulus'
readme.POWpre_selMNI_OCC = 'power around the max T value - Nvoxels  closest to the max modulated voxel, in prestimulus'
readme.POWpost_selMNI_OCC = 'power around the max T value - Nvoxels  closest to the max modulated voxel, in poststimulus'
readme.Occ_Cluster = '1 if significant occ cluster is present, 0 if not'
readme.Cluster = '1 if significant brain cluster is present, 0 if not'

readme.STAT = 'whole brain statistics in preselection of occipital areas '
readme.STAT_OCC = 'statistics'

save (filename, 'readme', 'Cluster', 'Occ_Cluster', 'MaxMNICoord', 'MaxMNI_Tstat', 'MaxisIn', 'MaxMNICoordOCC', 'MaxMNI_TstatOCC', 'POWpre_selMNI_OCC', 'POWpost_selMNI_OCC', 'MaxisInOCC', 'maxindOCC',...
      'POWpre_MaxMNI_OCC', 'POWpost_MaxMNI_OCC', 'POWpre_MaxMNI', 'POWpost_MaxMNI', 'STAT', 'STAT_OCC');

% d=sqrt(  (x1?x0).^2+(y1?y0).^2+(z1?z0).^2  );  % distance in 3D space

%%
exportToPPTX('saveandclose', PPTXname);



