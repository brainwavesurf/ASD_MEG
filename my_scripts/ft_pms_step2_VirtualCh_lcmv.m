% ft_step2_VirtualCh_lcmv
% cfg.method='lcmv' beamformer.subj
%
% Calculate spectral power in virtual sensors
% Calculate changes in the gamma range 

%For good solution this program uses the following:
% 1. Coregistration should be perfect: if source is in cerebellum, check the
% coregistration!
% 2. Use common covariance for all data  analysis [-1.0  1.2])
% 3. Use DFT filter on continious signal before epoching to remove 50 Hz
% artefact (0ptional). Not realy nesessary!

%%
clear
close all
screensize = get( groot, 'Screensize' );

% NB: add to path: 
fieldtripfolder = '/mnt/home/a_shishkina/fieldtrip/';
path(path, fieldtripfolder)
ft_defaults;
path('/mnt/home/a_shishkina/fieldtrip/external/mne/', path);

% This is an external matlab package used to save figures to PPTX
path('/mnt/home/a_shishkina/externals/pptx/', path);
path(path,'/mnt/home/a_shishkina/projects/asd_meg/0101/MEGdata/FT/'); %for fdr_bh.m

realDATAPATH = '/mnt/home/a_shishkina/projects/asd_meg/0101/MEGdata/';
DATAPATH  = '/mnt/home/a_shishkina/projects/asd_meg/0101/MEGdata/FT_beamf/';

tapsmofrq = 5; %25
gridres = 6; % 6 mm grid
lambda='5%';
alpha_betarange = [7 20];
%gammamax = [45 90]; % low gamma is often suppressed at fast velocity-V3, therefore, the power increase maximum will be defined for gamma >45Hz
Nvoxels = 25; % number of voxels [closest to the max] to average. We lloked at gamma spectra at this selection of voxels

filt = 'yes';
lpfreq=alpha_betarange(2)+5;
hpfreq=alpha_betarange(1)-5;

% colors for spectra
color{1}{1}=':k'; color{1}{2}='--k'; color{1}{3}='-k';
color{2}{1}=':b'; color{2}{2}='--b'; color{2}{3}='-b';
color{3}{1}=':g'; color{3}{2}='--g'; color{3}{3}='-g';


megfolder = strcat( '/meg', num2str(gridres), 'mm_nonlinwarp_tapsmofrq', num2str(tapsmofrq), '_lcmv', '/');
                              
%SUBJ = [ 'G001'; 'G001'; 'G002'; 'G002'; 'G003'; 'G003'; 'G005';'G006';'G007'; 'G007';  'G013'; 'G014'; 'G015'; 'G015'; 'G016'; 'G016'; 'G017'; 'G017'; 'G020'; 'G020'; 'G021'; 'G021'; 'G022'; 'G022'; 'G023'; 'G023'; 'G025'; 'G025'; 'G026'; 'G026'; 'G027'; 'G027'; 'G028'; 'G028'; 'G030'; 'G030'];
%DATE = [ 'f';    'l';  'f';    'l';       'f';    'l';    'l';   'f';   'f';   'l'; 'l'; 'l'; 'f'; 'l'; 'f'; 'l'; 'f'; 'l';  'f'; 'l'; 'f'; 'l'; 'f'; 'l'; 'f'; 'l';'f'; 'l'; 'f'; 'l'; 'f'; 'l'; 'f'; 'l'; 'f'; 'l'];
SUBJ = ['0101'];
DATE = ['oct-6'];
%% load atlas, template grid, common for all subjects

%  template MRI
templatefile = strcat (fieldtripfolder, '/external/spm8/templates/T1.nii'); 
template_mri = ft_read_mri(templatefile);
template_mri.coordsys = 'mni';

%  load atlas
atlas = ft_read_atlas( strcat (fieldtripfolder, '/template/atlas/aal/ROI_MNI_V4.nii') ); 
atlas = ft_convert_units(atlas,'cm'); % assure that atlas and template_grid are expressed in the %same units

% load template grid. Use the same template grid, that was used for the
% warped grid construction! Common for all subjects

load (strcat(DATAPATH, '0101_oct-6/mri_nonlinwarp_6mm_brthr0.5/0101_template_grid.mat'));
% load ( strcat(DATAPATH , subj, mrifolder, subj, '_template_grid.mat') );

cfg = [];
cfg.atlas      = atlas;
cfg.roi        = atlas.tissuelabel;  % here you can also specify a single label, i.e. single ROI
cfg.inputcoord = 'mni';
mask = ft_volumelookup(cfg, template_grid);

%chose only the 'brain tissue' 
template_grid.inside = false(template_grid.dim);
template_grid.inside(mask==1) = true;

% By doing this you get the label for every grid  :
cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';
atlas_template_grid = ft_sourceinterpolate(cfg, atlas, template_grid);
atlas_template_grid.tissue(isnan(atlas_template_grid.tissue)) =0;

ids      = find(atlas_template_grid.tissue);          % all interpolate regions
id       = atlas_template_grid.tissue(ids); %  all interpolate regions index
ROI      = atlas.tissuelabel(id);
occid1   = find(strncmpi(ROI,'Occipital',9));  %  indice
occid2   = find(strncmpi(ROI,'Calcarine',9));  %  indice
occid3   = find(strncmpi(ROI,'Cuneus',6));  %  indice
occid4   = find(strncmpi(ROI,'Lingual',7));  %  indice
occid5   = find(strncmpi(ROI,'Precuneus', 9));  
occid    = sort([occid1, occid2, occid3, occid4, occid5]);        
OCC      = ROI(occid);  % label

mask_occ = zeros(size(template_grid.pos,1), 1); mask_occ(ids(occid)) = 1;

%% across subjects

s=1;
size (SUBJ,1)
close all
subj = SUBJ(s,:);
date = DATE(s,:);

mkdir (strcat(DATAPATH, subj,'_', date, megfolder));
savemegto = strcat(DATAPATH, subj,'_', date, '/',megfolder);
mrifolder =strcat('/mri_nonlinwarp_', num2str(gridres), 'mm_', 'brthr',  '0.5/');

%% start PPTX report
exportToPPTX('new');
PPTXname  = strcat(savemegto, '/', subj,  '_lcmv_report');

%% subject specific Loads
% load leadfield / % source model    
load ( strcat(DATAPATH , subj,'_', date, '/', mrifolder, subj, '_grid_MNI_lf.mat') );
% head model: individ_hdm_vol
load ( strcat(DATAPATH ,subj,'_', date, '/', mrifolder, subj, '_individ_hdm_vol.mat' )); 
%     % DICS results
%     load (strcat(DATAPATH , subj, megfolder, subj, '_DISCbroad_source_maxT.mat'));
% load individual mri
load ( strcat(DATAPATH ,subj,'_', date, '/', mrifolder, subj, '_mri_orig_realigned.mat' )); 

%% read clean and delay corrected events exported from mne .py (exported with events_raw_resample_adults_4vel.py)

SUBJ2 = ['0102'];
subj2 = SUBJ2(s,:); 

raw_fiff_file = strcat(realDATAPATH, '0101_rings_ICA_raw.fif');
hdr = ft_read_header(raw_fiff_file);

ev_name=['/projects/asd_meg/0101/0101_clean_events.mat']
load(ev_name)
first= round(cast(hdr.orig.raw.first_samp, 'double'));
events(:,1) =  events(:,1)-first;

ev1 = find(events(:,3)==2);
ev2 = find(events(:,3)==4);
ev3 = find(events(:,3)==8);
EV = {ev1, ev2, ev3};

pre = -1.0* hdr.Fs ;
post = 1.2* hdr.Fs ;
trl=[];
for i=1:size (events,1)
    trl(i, 1)=(events(i,1)+pre) ; 
    trl(i, 2)=(events(i,1)+post) ; 
    trl(i, 3)= -1.0*hdr.Fs ; % offset
    trl(i, 4) = events(i,3); % stimulus_value;
end

% extract data epochs from the raw
cfg = [];
cfg.trl=trl;
cfg.channel     = 'meg';
cfg.dftfilter   = 'yes';
cfg.dftfreq     = [7 20];
cfg.demean = 'yes';
if strcmp(filt,'yes')  % filtering in the band-of-interest is recommended for lcmv
   cfg.lpfilter  = 'yes';
   cfg.lpfreq    = lpfreq;
   cfg.hpfilter  = 'yes';
   cfg.hpfreq    = hpfreq;
end
cfg.dataset = raw_fiff_file;
cfg    = ft_definetrial(cfg);
epochs = ft_preprocessing(cfg);

%% Plot average
cfg = [];
cfg.channel=epochs.label;
avg = ft_timelockanalysis(cfg,epochs);

hh=figure;
ttt = find(avg.time>-0.3 & avg.time<0.5);
plot (avg.time(ttt), avg.avg(:, ttt));
title ([subj, ': sensor average']);

exportToPPTX('addslide'); % slide with distributions
exportToPPTX('addpicture', hh, 'Position', [0.5,0.5,9,6]);
close (gcf);

%% Data covariance
cfg = [];
cfg.channel           = epochs.label;
cfg.covariance        = 'yes';
cfg.covariancewindow  =  'all'; % [-1.0  1.2];%'all'; %[-0.9 -0.1]; %
cfg.vartrllength      = 2;
total_avg = ft_timelockanalysis(cfg, epochs);
covariance.wind = cfg.covariancewindow;
covariance.epochs = 'all conditions';

%% perform source analysis
sourcemodel = grid_MNI_lf;
sourcemodel.inside = reshape(template_grid.inside, [1, template_grid.dim(1)*template_grid.dim(2)*template_grid.dim(3)]);

cfg=[];
cfg.method='lcmv';
cfg.grid = sourcemodel; %template_grid;%grid_MNI_lf;
cfg.headmodel=individ_hdm_vol;
cfg.lcmv.keepfilter='yes'; % keep filters in the output, which are later multiplied with the data
cfg.lcmv.fixedori='yes'; % id 'yes' consider only the dominant orientation
cfg.lcmv.lambda=lambda;        
%cfg.lcmv.projectmom = 'yes';
cfg.reducerank = 2;
cfg.normalize = 'yes'; %leadfield is computed on the fly. In this case you can cormalise it
source_total=ft_sourceanalysis(cfg, total_avg);
source_total.pos = template_grid.pos;

%% 
for con = 1:3 % for conditions

% select epochs
cfg = [];
cfg.trials =  EV{con};        
EPO = ft_selectdata(cfg, epochs);

cfg = [];
cfg.channel=EPO.label;
avg = ft_timelockanalysis(cfg,EPO);
AVG{con}=avg;

% subtract evoked! Important for Moscow data with 60Hz projector.
for ttt=1:size (EPO.trial,2)
    EPO.trial{ttt}=EPO.trial{ttt}-avg.avg;
end

% select pre and poststim   
cfg=[];
cfg.latency     = [-0.9 0]; EPOCH.prestimulus=cfg.latency;
[data2_pre] = ft_selectdata(cfg, EPO);

cfg.latency     = [0.3 1.2]; EPOCH.poststimulus=cfg.latency;
[data2_post] = ft_selectdata(cfg, EPO);

%% This is for cfg.lcmv.fixedori='yes'; 
% Multiply filters with the data and organize into FieldTrip sensable data structure
spatialfilter_pre=cat(1,source_total.avg.filter{:});  % all voxels inside
spatialfilter_post=cat(1,source_total.avg.filter{:});  % all voxels inside

virtsens_pre=[]; virtsens_post=[];
for i=1:length(data2_pre.trial)
    virtsens_pre.trial{i}=spatialfilter_pre*data2_pre.trial{i};
    virtsens_post.trial{i}=spatialfilter_post*data2_post.trial{i};
end
virtsens_pre.time=data2_pre.time;
virtsens_pre.fsample=data2_pre.fsample;
virtsens_pre.label= cellstr(string(find(sourcemodel.inside))); %(occid(isnotempt)))'
virtsens_pre.pos = template_grid.pos(find(sourcemodel.inside),:);   
virtsens_pre.grad =data2_pre.grad;

virtsens_post.time=data2_post.time;
virtsens_post.fsample=data2_post.fsample;
virtsens_post.label= cellstr(string(find(sourcemodel.inside))); %(occid(isnotempt)))'
virtsens_post.pos = template_grid.pos(find(sourcemodel.inside),:);   
virtsens_post.grad =data2_post.grad;

%% do spectral analysis
cfg = [];
cfg.method    = 'mtmfft';
cfg.output    = 'pow';
cfg.foilim    = [5 alpha_betarange(2)];
cfg.pad = 'nextpow2';
cfg.tapsmofrq = tapsmofrq;
cfg.keeptrials = 'yes';
freq1          = ft_freqanalysis(cfg, virtsens_pre); % keeptrial
freq2          = ft_freqanalysis(cfg, virtsens_post);
cfg.keeptrials   = 'no';
fd_pre            = ft_freqdescriptives(cfg, freq1); % average
fd_post            = ft_freqdescriptives(cfg, freq2);

% %     %% Spectrum in the maximally modulated voxel plus surround
% %     % load saved: maxindOCC from 
% %     loadfilename = strcat(DATAPATH, subj, '/',  'meg6mm_linwarp_tapsmofrq5_individWF_maxchange_static/', subj, '_DISCbroad_source_maxT.mat');
% %     load (loadfilename, 'maxindOCC'); % max probability voxel by DICS
% %     
% % % % maxindOCC   - index with max probability: max modulated voxel
% %     inside_id = find(template_grid.inside);
% %     max_id =find (inside_id == maxindOCC{con})
% %     inside_max_pos = virtsens_post.pos(max_id ,:);   
% %     
% %     % find 6 mm distance
% %     x0 = inside_max_pos(1); y0 = inside_max_pos(2); z0 = inside_max_pos(3);
% %     d= sqrt(  (virtsens_post.pos(:,1)-x0).^2+(virtsens_post.pos(:,2)-y0).^2+(virtsens_post.pos(:,3)-z0).^2 );
% %     A = [d, inside_id, [1:length(d)]' ];
% %     B=sortrows(A,1);
% %     n=1:16;
% %     take_voxels_inside = B(n,3);
% %     virtsens_post.pos(take_voxels_inside,:);
% %    
% %     AA=[d, inside_id,  virtsens_post.pos];
% %     
%%
f1=gammamax(1); f2=gammamax(2);
gamma_ind = find(fd_pre.freq>=floor(f1) & fd_pre.freq<=ceil(f2) ); % gamma frequencues
gamma_wholerange_ind = find(fd_pre.freq>=floor(gammarange(1)) & fd_pre.freq<=ceil(gammarange(2)) ); % gamma frequencues
inside_id = find(template_grid.inside);
DIFFgamma =  (mean( fd_post.powspctrm(:, gamma_ind),2)-mean(fd_pre.powspctrm(:, gamma_ind),2))./mean(fd_pre.powspctrm(:, gamma_ind) ,2); % average diff in gamma band

Spectrum_pre{con} = fd_pre;
Spectrum_post{con} = fd_post;

% plot DIFFgamma
SourseDIFF = source_total;
SourseDIFF.pow = DIFFgamma;
SourseDIFF = rmfield(SourseDIFF,'avg')
SourseDIFF.pow = NaN(size(SourseDIFF.inside))';
SourseDIFF.pow(find(SourseDIFF.inside)) = DIFFgamma';
SourseDIFF.coordsys = 'spm'

% interpolate to template brain 
[SourseDIFF] = ft_convert_units(SourseDIFF,'cm')
[template] = ft_convert_units(template_mri,'cm')
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
SourseDIFF_Int  = ft_sourceinterpolate(cfg, SourseDIFF, template);
SourseDIFF_Int.coordsys= 'mni';

cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'pow';
cfg.maskparameter = cfg.funparameter;
%cfg.funcolorlim   = [0.0 0.005];
%cfg.location = [0, -10, 0];
%cfg.opacitylim    = 'zeromax'; 
cfg.opacitymap    = 'rampup'; 
cfg.atlas=atlas;
%cfg.roi = {'Occipital_Sup_R', 'Occipital_Sup_L',  'Occipital_Inf_R',  'Occipital_Inf_L', 'Cuneus_R', 'Cuneus_L', 'Precuneus_R' , 'Precuneus_L', 'Occipital_Mid_R' , 'Occipital_Mid_L', 'Calcarine_R', 'Calcarine_L', 'Lingual_R' , 'Lingual_L'      };
ft_sourceplot(cfg, SourseDIFF_Int );

% save distribution of the ratio:
exportToPPTX('addslide'); % slide with distributions
exportToPPTX('addpicture', gcf, 'Position', [0.3,1,4,3]);
exportToPPTX('addtext', strcat(subj, ': lcmv beamformer, total and occipital, condition V', num2str(con), ': (post-pre)/pre in the gamma range: ', num2str(gammamax(1)), '-', num2str(gammamax(2)),'Hz, filter=', filt  ), 'FontSize', 12, 'Position', [0.5, 0.5, 10,1]);
close(gcf);

% Plot with occ mask
SourseDIFF.inside = zeros(size(SourseDIFF.inside));
SourseDIFF.inside(ids(occid))  = 1;
[SourseDIFF] = ft_convert_units(SourseDIFF,'cm')
[template] = ft_convert_units(template_mri,'cm')
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
SourseDIFF_Int  = ft_sourceinterpolate(cfg, SourseDIFF, template);
SourseDIFF_Int.coordsys= 'mni';
cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'pow';
cfg.maskparameter = cfg.funparameter;
%cfg.funcolorlim   = [0.0 0.005];
%cfg.location = [0, -10, 0];
%cfg.opacitylim    = 'zeromax'; 
cfg.opacitymap    = 'rampup'; 
cfg.atlas=atlas;
ft_sourceplot(cfg, SourseDIFF_Int );
exportToPPTX('addpicture', gcf, 'Position', [4.5,1,4,3]);
close(gcf);

%% Save positions and probability for maximally modulated voxel, whole brain
[DIFFmax,x] = max(DIFFgamma);    

A = [DIFFgamma, inside_id, [1:length(inside_id)]', virtsens_post.pos ];
B=flip(sortrows(A,1),1); 
max_voxel_inside = B(1,3); %index inside
max_voxel_pos = virtsens_post.pos(max_voxel_inside,:);
MAX_VOXEL_NUMBER{con} = max_voxel_inside;

%%%% for plot 
n=1:Nvoxels;
max_voxels_inside = B(n,3);
max_voxels_pos = virtsens_post.pos(max_voxels_inside,:);
MAX_ROI{con}= ROI(B(n,3));
MAX_coord{con}= source_total.pos(B(n,2),:);
%%%%%

% probability of in the max modulated voxel
pre =  squeeze(mean(freq1.powspctrm(:, max_voxel_inside, gamma_wholerange_ind),2)); % tr, voxel_inside, freq; / average over max_voxels_inside and gamma_ind
post = squeeze(mean(freq2.powspctrm(:, max_voxel_inside, gamma_wholerange_ind),2));
for j=1: size(pre,2) % for all freq
    [h, p] = ttest(squeeze(pre(:,j)),squeeze(post(:,j)))
    PValues(j)=p;
end
[FDR05max{con}.h, FDR05max{con}.crit_p, FDR05max{con}.adj_p]=fdr_bh(PValues ,0.0001 ,'pdep','yes');
[FDR0001max{con}.h, FDR0001max{con}.crit_p, FDR0001max{con}.adj_p]=fdr_bh(PValues ,0.0001 ,'pdep','yes');
maxvoxel_pre{con} = fd_pre.powspctrm(max_voxel_inside,:);
maxvoxel_post{con} = fd_post.powspctrm(max_voxel_inside,:);
MAX_ROI{con}= ROI(B(1,3));
MAX_coord{con}= source_total.pos(B(1,2),:);

% p for averaged gammamax range
pre =  squeeze(mean(mean(freq1.powspctrm(:, max_voxel_inside, gamma_ind),3),2)); % tr, voxel_inside, freq; / average over max_voxels_inside and gamma_ind
post = squeeze(mean(mean(freq2.powspctrm(:, max_voxel_inside, gamma_ind),3),2));
[Pgamma_max{con}, h, statis] = signrank(pre,post); % P in the gamma range in the max voxel

%%  Save positions and probability for maximally modulated voxel, occipital region
% ids, occid
Aocc = [DIFFgamma(occid), inside_id(occid), occid', virtsens_post.pos(occid,:) ];
Bocc=flip(sortrows(Aocc,1),1);
max_voxel_inside_occ = Bocc(1,3);
max_voxel_pos_occ = virtsens_post.pos(max_voxel_inside_occ,:);
MAX_VOXEL_OCC_NUMBER{con} = max_voxel_inside_occ;

%%%% for plot 
n=1:Nvoxels;
max_voxels_inside_occ = Bocc(n,3);
max_voxels_pos_occ = virtsens_post.pos(max_voxels_inside_occ,:);
MAX_ROIocc{con}= ROI(Bocc(n,3));
MAX_coord_occ{con}= source_total.pos(Bocc(n,2),:);
%%%%%

% probability in the max modulated occ voxels
pre_occ =  squeeze(freq1.powspctrm(:, max_voxel_inside_occ, gamma_wholerange_ind) ); % tr, voxel_inside, freq; / average over max_voxels_inside and gamma_ind
post_occ = squeeze(freq2.powspctrm(:, max_voxel_inside_occ, gamma_wholerange_ind) );
for j=1: size(pre_occ,2) % for all freq
    [h, p] = ttest(squeeze(pre_occ(:,j)),squeeze(post_occ(:,j)));
    PValues(j)=p;
end
[FDR05occmax{con}.h, FDR05occmax{con}.crit_p, FDR05occmax{con}.adj_p]=fdr_bh(PValues ,0.0001 ,'pdep','yes');
[FDR0001occmax{con}.h, FDR0001occmax{con}.crit_p,  FDR0001occmax{con}.adj_p]=fdr_bh(PValues ,0.0001 ,'pdep','yes');
maxvoxel_pre_occ{con} = fd_pre.powspctrm(max_voxel_inside_occ,:);
maxvoxel_post_occ{con} = fd_post.powspctrm(max_voxel_inside_occ,:);
MAX_ROIocc{con}= ROI(Bocc(1,3));
MAX_coord_occ{con}= source_total.pos(Bocc(1,2),:);

% p for averaged gammamax range
pre =  squeeze(mean(freq1.powspctrm(:, max_voxel_inside_occ, gamma_ind),3)); % tr, voxel_inside, freq; / average over max_voxels_inside and gamma_ind
post = squeeze(mean(freq2.powspctrm(:, max_voxel_inside_occ, gamma_ind),3));
[Pgamma_occmax{con}, h, statis] = signrank(pre,post); % P in the gammamax range in the occipital max voxel, for averaged gammamax range

PREoccmax{con}  = squeeze(fd_pre.powspctrm(max_voxel_inside_occ, :));
POSToccmax{con} = squeeze(fd_post.powspctrm(max_voxel_inside_occ, :));
DIFFoccmax{con} = (POSToccmax{con}-PREoccmax{con})./PREoccmax{con};

% wF and power for DIFFoccmax
gammaind = find(fd_post.freq>=gammarange(1) &  fd_post.freq<=gammarange(2));
[maxPow{con},  ind] = max(DIFFoccmax{con}(gammaind));
maxF{con} = fd_post.freq(gammaind(ind));
inds = find(DIFFoccmax{con}(gammaind) >=2/3*maxPow{con});
% if frequencies are too far from the maximum (more than 20 Hz) we exclude them!
excl = find( (fd_post.freq(gammaind(inds))-maxF{con}) >20);
inds(excl) =[];
binsfreqsmax{con} = fd_post.freq(gammaind(inds));

maxWF{con} = sum(fd_post.freq(gammaind(inds)).*DIFFoccmax{con}(gammaind(inds) ))/sum(DIFFoccmax{con}(gammaind(inds) ));
maxWPow{con} = mean(DIFFoccmax{con}(gammaind(inds)) );

%FDR05occmax, FDR0001occmax,    PREoccmax, POSToccmax, DIFFoccmax,     Pgamma_max,      maxPow, maxF, maxWF, maxWPow

%% average power in voxels closest to maxmodulated OCC voxel: Nvoxels averaged
% calculate distance from the max voxel  MAX_voxel_inside_occ  DIFFgamma
% fd_post fd_pre  virtsens_post.pos
x1=virtsens_post.pos(:,1); y1=virtsens_post.pos(:,2);  z1=virtsens_post.pos(:,3); 
x0=virtsens_post.pos(max_voxel_inside_occ,1); y0=virtsens_post.pos(max_voxel_inside_occ,2);  z0=virtsens_post.pos(max_voxel_inside_occ,3); 
distance = sqrt(  (x1-x0).^2 +(y1-y0).^2 +(z1-z0).^2 ); 
distance = [(1:length(distance))' , distance];
Cocc=sortrows(distance,2);
Nvoxels_id =Cocc(1:Nvoxels,:); Nvoxels_id = Nvoxels_id(:,1);
NVOXELS_ID{con} = Nvoxels_id;

PREoccsel{con} = mean(fd_pre.powspctrm(Nvoxels_id, :),1);
POSToccsel{con} = mean(fd_post.powspctrm(Nvoxels_id, :),1);
DIFFoccsel{con} = (POSToccsel{con}-PREoccsel{con})./PREoccsel{con};

% wF and power for DIFFoccsel
% probability in the max selection occ voxels
pre_occ =  squeeze(mean(freq1.powspctrm(:,  Nvoxels_id, gamma_wholerange_ind),2)); % tr, voxel_inside, freq; / average over max_voxels_inside and gamma_ind
post_occ = squeeze(mean(freq2.powspctrm(:,  Nvoxels_id, gamma_wholerange_ind),2));
for j=1: size(pre_occ,2) % for all freq
    [h, p] = ttest(squeeze(pre_occ(:,j)),squeeze(post_occ(:,j)));
     PValues(j)=p;
end
[FDR05occsel{con}.h, FDR05occsel{con}.crit_p, FDR05occsel{con}.adj_p]=fdr_bh(PValues ,0.0001 ,'pdep','yes');
[FDR0001occsel{con}.h, FDR0001occsel{con}.crit_p, FDR0001occsel{con}.adj_p]=fdr_bh(PValues ,0.0001 ,'pdep','yes');
selvoxel_pre_occ{con}  = mean(fd_pre.powspctrm(Nvoxels_id,:),1);
selvoxel_post_occ{con} = mean(fd_post.powspctrm(Nvoxels_id,:),1);

% p for averaged gammamax range
pre =  squeeze(mean(mean(freq1.powspctrm(:, Nvoxels_id, gamma_ind),3),2)); % tr, voxel_inside, freq; / average over max_voxels_inside and gamma_ind
post = squeeze(mean(mean(freq2.powspctrm(:, Nvoxels_id, gamma_ind),3),2));
[Pgammarange_sel{con}, h, statis] = signrank(pre,post); % P in the gamma range in the occipital Nvoxels

PREoccsel{con}  = mean(fd_pre.powspctrm(Nvoxels_id, :),1);
POSToccsel{con} = mean(fd_post.powspctrm(Nvoxels_id, :),1);
DIFFoccsel{con} = (POSToccsel{con}-PREoccsel{con})./PREoccsel{con};

% wF and power for DIFFoccsel
[selPow{con},  ind] = max(DIFFoccsel{con}(gammaind));
selF{con} = fd_post.freq(gammaind(ind));
inds = find(DIFFoccsel{con}(gammaind) >=2/3*selPow{con});
% if frequencies are too far from the maximum (more than 20 Hz) we exclude them!
excl = find( (fd_post.freq(gammaind(inds))-selF{con}) >20);
inds(excl) =[];
binsfreqssel{con} = fd_post.freq(gammaind(inds));
INDS{con}=inds; % gamma band indexes to remember

selWF{con} = sum(fd_post.freq(gammaind(inds)).*DIFFoccsel{con}(gammaind(inds) ))/sum(DIFFoccsel{con}(gammaind(inds) ));
selWPow{con} = mean(DIFFoccsel{con}(gammaind(inds)) );

% FDR05occsel, FDR0001occsel,    PREoccsel, POSToccsel, DIFFoccsel,     Pgamma_sel,     selPow, selF, selWF, selWPow

%% Save Source info:
SORTED{con}=B;
SORTEDocc{con}=Bocc;

%   %% In each condition plot spectrum for max modulated voxels in occipital selection
%    POST=mean(fd_post.powspctrm(max_voxels_inside_occ,:),1);
%    PRE=mean(fd_pre.powspctrm(max_voxels_inside_occ,:),1);
%    DIFF{con} = (POST-PRE)./PRE;
% 

%%  Plot spectra for 10 the maximally modutated voxels, whole brain 
N= 10;
pos_fig1 = [10 10 screensize(3)/10*9 screensize(4)/10*9];    
hh=figure('Position',pos_fig1);
for i=1:N
    i
    POST=mean(fd_post.powspctrm(max_voxels_inside(i),:),1);
    PRE=mean(fd_pre.powspctrm(max_voxels_inside(i),:),1);
    diff = (POST-PRE)./PRE;
    subplot(2,N/2,i);
    p1=plot (fd_pre.freq(gamma_wholerange_ind), PRE(gamma_wholerange_ind), 'b');  hold on
    p2=plot (fd_pre.freq(gamma_wholerange_ind), POST(gamma_wholerange_ind), 'r'); 
    xlabel('Frequency, Hz')
    ylabel('Power')
    yyaxis right
    ylabel('(post-pre)/pre')
    p3=plot(fd_pre.freq(gamma_wholerange_ind), diff(gamma_wholerange_ind), '--k');
    hold off
    title (  strcat ( ROI(B(i,3)), ', xyz:', num2str(B(i,4:6), '%4.1f ') ) ) ;        
end

lgd=legend([p1 p2 p3], {'pre pow', 'post pow', '(post-pre)/pre'});
lgd.FontSize = 8;
suptitle(strcat( subj, ': pre- and post-stimulus power and (post-pre)/pre ratio in gamma range in condition..', num2str(con), ', 10 max voxels in the whole brain') )
exportToPPTX('addslide'); % slide with distributions
exportToPPTX('addpicture', hh, 'Position', [0.5,1,9,6]);
exportToPPTX('addtext', strcat('WHOLE BRAIN: ', subj, ', Condition:', num2str(con),  ', Filter:', filt, ', gamma range:', num2str(gammarange(1)), '-',...
              num2str(gammarange(2)), 'Hz ', ', cfg.lcmv.lambda=', lambda, ', tapsmofrq=', num2str(tapsmofrq),  ', Ranksun p for the mean power in 45-110Hz:', num2str(Pgamma_max{con})   ),  'FontSize', 14, 'Position', [0.1, 0.1, 10,1]);
close(gcf);

%%  Plot spectra for 10 the maximally modutated voxels, occipital selection
pos_fig1 = [10 10 screensize(3)/10*9 screensize(4)/10*9];    
hh=figure('Position',pos_fig1);
for i=1:N
    i
    POST=mean(fd_post.powspctrm(max_voxels_inside_occ(i),:),1);
    PRE=mean(fd_pre.powspctrm(max_voxels_inside_occ(i),:),1);
    diff = (POST-PRE)./PRE;
    subplot(2,N/2,i);
    p1=plot (fd_pre.freq(gamma_wholerange_ind), PRE(gamma_wholerange_ind), 'b');  hold on
    p2=plot (fd_pre.freq(gamma_wholerange_ind), POST(gamma_wholerange_ind), 'r'); 
    xlabel('Frequency, Hz')
    ylabel('Power')
    yyaxis right
    ylabel('(post-pre)/pre')
    p3=plot(fd_pre.freq(gamma_wholerange_ind), diff(gamma_wholerange_ind), '--k');
    hold off
    title (  strcat (ROI(Bocc(i,3)), ', xyz:', num2str(Bocc(i,4:6), '%4.1f ') ) ) ;        
end
lgd=legend([p1 p2 p3], {'pre', 'post', '(post-pre)/pre'}, 'FontSize', 16);
lgd.FontSize = 8;
suptitle(strcat('Subject:', subj, ', pre- and post-stimulus power and (post-pre)/pre ratio in gamma range in condition..', num2str(con), ', 10 maximal voxels in occipital selection') )

exportToPPTX('addslide'); % slide with distributions
exportToPPTX('addpicture', hh, 'Position', [0.5,1,9,6]);
exportToPPTX('addtext', strcat('OCCIPITAL SELECTION:   Subject:', subj, ', Condition:', num2str(con),  ', Filter:', filt, ', gamma range:', num2str(gammarange(1)), '-',...
                 num2str(gammarange(2)), 'Hz ', ', cfg.lcmv.lambda=', lambda, ', tapsmofrq=', num2str(tapsmofrq), ', Minimal FDR-p-level for max-modulated-voxel in gamma range =', num2str(min(FDR05occmax{con}.adj_p)), ', Ranksun p for the mean power in 45-110Hz:', num2str(Pgamma_occmax{con})  ), 'FontSize', 14, 'Position', [0.1, 0.1, 10,1]);
close(gcf);

%%
FREQ1powspctrm{con} = squeeze(mean(freq1.powspctrm(:,:,gamma_ind),3)); %keeptrial
FREQ2powspctrm{con} = squeeze(mean(freq2.powspctrm(:,:,gamma_ind),3));

end % end for conditons

% % %     %% power in 25 max voxels based on 100% contrast condition, for contrast 50%
% % %     % load the voxel selection from the subject's 100% contrast experiment
% % %     LC1 = load (strcat('/Volumes/MEGdisk/Moscow_4velocities/Beamformer/', subj, '/meg6mm_linwarp_tapsmofrq5_lcmv_static/', subj, '_static_Lcmv_source_spectra.mat'),  'selWPow'); 
% % %     LC2 = load (strcat('/Volumes/MEGdisk/Moscow_4velocities/Beamformer/', subj, '/meg6mm_linwarp_tapsmofrq5_lcmv_static/', subj, '_static_Lcmv_source_spectra.mat'), 'NVOXELS_ID');
% % %     [a,b] = max(cell2mat(LC1.selWPow)); % b is the max condition
% % %        for con=1:4  
% % %            pre =  squeeze(mean(FREQ1powspctrm{con}(:, NVOXELS_ID{b}),2)); % tr, voxel_inside, freq; / average over max_voxels_inside and gamma_ind
% % %            post = squeeze(mean(FREQ2powspctrm{con}(:, NVOXELS_ID{b}),2));
% % %            [Pgammarange_sel_MAXselection{con}, h, statis] = signrank(pre,post); % P in the gamma range in the occipital Nvoxels
% % % 
% % %            PREoccsel_MAXselection{con}  = mean(Spectrum_pre{con}.powspctrm(NVOXELS_ID{b}, :),1);
% % %            POSToccsel_MAXselection{con} = mean(Spectrum_post{con}.powspctrm(NVOXELS_ID{b}, :),1);
% % %            DIFFoccsel_MAXselection{con} = (POSToccsel_MAXselection{con}-PREoccsel_MAXselection{con})./PREoccsel_MAXselection{con};
% % % 
% % %            % wFselMAX and power for DIFFoccsel_MAXselection
% % %            [selMAXPow{con},  ind] = max(DIFFoccsel{con}(gammaind));
% % %            selMAXF{con} = fd_post.freq(gammaind(ind));
% % %            inds = find(DIFFoccsel{con}(gammaind) >=2/3*selMAXPow{con});
% % %             % if frequencies are too far from the maximum (more than 20 Hz) we exclude them!
% % %            excl = find( (fd_post.freq(gammaind(inds))-selMAXF{con}) >20);
% % %            inds(excl) =[];
% % %            binsfreqsselMAX{con} = fd_post.freq(gammaind(inds));
% % % 
% % %            selWMAXF{con} = sum(fd_post.freq(gammaind(inds)).*DIFFoccsel_MAXselection{con}(gammaind(inds) ))/sum(DIFFoccsel_MAXselection{con}(gammaind(inds) ));
% % %            selWMAXPow{con} = mean(DIFFoccsel_MAXselection{con}(gammaind(inds)) );
% % %        end  
% % % 
% % %    end
%% power in 25 max voxels based on the MAX condition:  
[a,b] = max(cell2mat(selWPow)); % b is the max condition
 for con=1:3  
     pre =  squeeze(mean(FREQ1powspctrm{con}(:, NVOXELS_ID{b}),2)); % tr, voxel_inside, freq; / average over max_voxels_inside and gamma_ind
     post = squeeze(mean(FREQ2powspctrm{con}(:, NVOXELS_ID{b}),2));
     [Pgammarange_sel_MAXselection{con}, h, statis] = signrank(pre,post); % P in the gamma range in the occipital Nvoxels

     PREoccsel_MAXselection{con}  = mean(Spectrum_pre{con}.powspctrm(NVOXELS_ID{b}, :),1);
     POSToccsel_MAXselection{con} = mean(Spectrum_post{con}.powspctrm(NVOXELS_ID{b}, :),1);
     DIFFoccsel_MAXselection{con} = (POSToccsel_MAXselection{con}-PREoccsel_MAXselection{con})./PREoccsel_MAXselection{con};

     % wFselMAX and power for DIFFoccsel_MAXselection
     [selMAXPow{con},  ind] = max(DIFFoccsel{con}(gammaind));
     selMAXF{con} = fd_post.freq(gammaind(ind));
     inds = find(DIFFoccsel{con}(gammaind) >=2/3*selMAXPow{con});
     % if frequencies are too far from the maximum (more than 20 Hz) we exclude them!
     excl = find( (fd_post.freq(gammaind(inds))-selMAXF{con}) >20);
     inds(excl) =[];
     binsfreqsselMAX{con} = fd_post.freq(gammaind(inds));

     selWMAXF{con} = sum(fd_post.freq(gammaind(inds)).*DIFFoccsel_MAXselection{con}(gammaind(inds) ))/sum(DIFFoccsel_MAXselection{con}(gammaind(inds) ));
     selWMAXPow{con} = mean(DIFFoccsel_MAXselection{con}(gammaind(inds)) );

 end 

%% save figure with spectra in max selection
figure;
plot(fd_pre.freq, PREoccsel{1}, color{1}{1}); hold on; plot(fd_post.freq,POSToccsel{1},color{1}{2});   
plot(fd_pre.freq, PREoccsel{2}, color{2}{1});          plot(fd_post.freq,POSToccsel{2},color{2}{2});   
plot(fd_pre.freq, PREoccsel{3}, color{3}{1});          plot(fd_post.freq,POSToccsel{3},color{3}{2});   



xlabel('Frequency, Hz');
ylabel('Power')
yyaxis right
ylabel('Difference: (pre-post)/pre')
plot(fd_post.freq, DIFFoccsel{1}, color{1}{3}); 
plot(fd_post.freq, DIFFoccsel{2}, color{2}{3}); 
plot(fd_post.freq, DIFFoccsel{3}, color{3}{3}); 

hold off;

lgd=legend( {'SlowPre', 'SlowPost', 'MediumPre', 'MediumPost', 'FastPre', 'FastPost', 'MediumDiff','FastDiff' });
lgd.FontSize = 9;

title (strcat('Ave of..', num2str(Nvoxels), '..vox closest to occ max'), 'FontSize', 12)
exportToPPTX('addslide'); % slide with distributions
exportToPPTX('addpicture', gcf, 'Position', [1,1,7,4.5]);
exportToPPTX('addtext', strcat('Subj:', subj, ', spectral changes in selection of...', num2str(Nvoxels), ' voxels, LCMV beamformer' ), 'FontSize', 12, 'Position', [0.5, 0.5, 10,1]);
close(gcf);

%% save figure with spectra in the max condition selection
figure;
plot(fd_pre.freq, PREoccsel_MAXselection{1}, color{1}{1}); hold on; plot(fd_post.freq,POSToccsel_MAXselection{1},color{1}{2}); 
plot(fd_pre.freq, PREoccsel_MAXselection{2}, color{2}{1});          plot(fd_post.freq,POSToccsel_MAXselection{2},color{2}{2}); 
plot(fd_pre.freq, PREoccsel_MAXselection{3}, color{3}{1});          plot(fd_post.freq,POSToccsel_MAXselection{3},color{3}{2}); 

xlabel('Frequency, Hz');
ylabel('Power')
yyaxis right
ylabel('Difference: (pre-post)/pre')
plot(fd_post.freq, DIFFoccsel_MAXselection{1}, color{1}{3}); 
plot(fd_post.freq, DIFFoccsel_MAXselection{2}, color{2}{3}); 
plot(fd_post.freq, DIFFoccsel_MAXselection{3}, color{3}{3}); 
hold off;

lgd=legend( {'SlowPre', 'SlowPost', 'MediumPre', 'MediumPost', 'FastPre', 'FastPost', 'SlowDiff','MediumDiff','FastDiff' });
lgd.FontSize = 9;
title (strcat('Ave..', num2str(Nvoxels), '..vox closest to occ max in MAX cond.'), 'FontSize', 12)
exportToPPTX('addslide'); % slide with distributions
exportToPPTX('addpicture', gcf, 'Position', [1,1,7,4.5]);
exportToPPTX('addtext', strcat('Subj:', subj, ', spectral changes in selection of...', num2str(Nvoxels), ' voxels, LCMV beamformer' ), 'FontSize', 12, 'Position', [0.5, 0.5, 10,1]);
close(gcf);

%% save figure with spectra in max 
figure;
plot(fd_pre.freq, PREoccmax{1}, color{1}{1}); hold on; plot(fd_post.freq,POSToccmax{1},color{1}{2}); 
plot(fd_pre.freq, PREoccmax{2}, color{2}{1});          plot(fd_post.freq,POSToccmax{2},color{2}{2}); 
plot(fd_pre.freq, PREoccmax{3}, color{3}{1});          plot(fd_post.freq,POSToccmax{3},color{3}{2}); 

xlabel('Frequency, Hz');
ylabel('Power')
yyaxis right
ylabel('Difference: (pre-post)/pre')
plot(fd_post.freq, DIFFoccmax{1}, color{1}{3});
plot(fd_post.freq, DIFFoccmax{2}, color{2}{3}); 
plot(fd_post.freq, DIFFoccmax{3}, color{3}{3});
hold off;

lgd=legend( {'SlowPre', 'SlowPost', 'MediumPre', 'MediumPost', 'FastPre', 'FastPost', 'SlowDiff','MediumDiff','FastDiff' });
lgd.FontSize = 9;
title ('1 Max occipital voxel ', 'FontSize', 12)
exportToPPTX('addslide'); % slide with distributions
exportToPPTX('addpicture', gcf, 'Position', [1,1,7,4.5]);
exportToPPTX('addtext', strcat('Subj:', subj, ', spectral changes in selection of...', num2str(Nvoxels), ' voxels, LCMV beamformer' ), 'FontSize', 12, 'Position', [0.5, 0.5, 10,1]);
close(gcf);
%% Save spectra and diff
%filename = strcat(savemegto, subj, '_Lcmv_source_spectra.mat');
readme = []
readme.covariance = 'covatiance window and epochs info'
readme.SORTED = 'sorted by diff=post-pre/pre:  [diff, voxel, voxel_inside, X-position, Y-position, Z-position  ], positions in MNI coordinates'
readme.info = strcat('Subject:', subj, ', Condition:', num2str(con),  ', Filter:', filt, ', gamma range:', num2str(f1), '-',...
                 num2str(f2), 'Hz ', ', cfg.lcmv.lambda=', lambda, ', tapsmofrq=', num2str(tapsmofrq)  );
readme.freq = 'frequencies';
readme.FDR = 'false discovery rate corrected p at the [average of] maximally modulated voxels, with 05 and 0001 thresholds';

readme.PREocc  = 'whole power spectrum in the max modulated voxels or selection (of Nxovels), prestimulus';
readme.POSTocc = 'whole power spectrum in the max modulated voxels or selection (of Nxovels), poststimulus';
readme.DIFFocc = 'whole power spectrum in the max modulated voxels or selection (of Nxovels), (post-pre)/pre difference';

readme.Pgamma_max       = 'P in the gammamax range in the           max voxel, for averaged gammamax range: Ranksum p for the gammamax 45-90 Hz range';
readme.Pgamma_occmax    = 'P in the gammamax range in the occipital max voxel, for averaged gammamax range: Ranksum p for the gammamax 45-90 Hz range';
readme.Pgammarange_sel  = 'probability of increase  in the averaged gamma range at max-occ selection of 25 voxels: Ranksum p for the gammamax 45-90 Hz range';
readme.Pow   = 'max power for the max voxel/selection of Nvoxels, calculated from normalized power of type (POSToccmax{con}-PREoccmax{con})./PREoccmax{con}'  
readme.F     = 'freq of the max for the max voxel/selection of Nvoxels' 
readme.WPow  = 'weighted power for the max voxel/selection of Nvoxels, calculated from normalized power of type (POSToccmax{con}-PREoccmax{con})./PREoccmax{con}'
readme.WF    = 'weighted freq for the max voxel/selection of Nvoxels'

readme.selWMAXF = 'the same as selWF, but the selection of N voxels is based on MAX condition'
readme.selWMAXPow = 'the same as selWPow, but the selection of N voxels is based on MAX condition'
readme.Pgammarange_sel_MAXselection = 'probability of increase  in the averaged gamma range at max-occ selection of 25 voxelsvased on MAX: Ranksum p for the gammamax 45-90 Hz range';

readme.MAX_ROI = 'region of the maximally modulated voxel: in the whole brain or in occipital selection - occ'
readme.MAX_coord = 'coordinates of N maximal voxels in whole brain/occipital'

readme.Spectrum_pre_and_post = 'full brain spectra, pre and post-stimulus';
readme.DIFF = '(post-pre)/pre at each frequency in selection of Nvoxels max occ voxels'
readme.binsfreqssel = 'bins used wor selWF or maxWF'

readme.NVOXELS_ID = 'id nimbers of 25 voxels around the max voxel';

readme.MAX_VOXEL_OCC_NUMBER ='number of maximally modulated voxel in occipital selection';

readme.PREvisual  =  'mean normalized power in visual areas in prestimulus';
readme.POSTvisual =  'mean normalized power in visual areas in poststimulus';
readme.PREvisual  =  'mean normalized power in calcarines in prestimulus';
readme.POSTvisual =  'mean normalized power in calcarines in poststimulus';
readme.EPOCH = 'epoch size';
Freq = freq1.freq; % gamma range frequencies

    % FDR05occsel, FDR0001occsel,    PREoccsel, POSToccsel, DIFFoccsel,     Pgamma_sel,     selPow, selF, selWF, selWPow
    % FDR05occmax, FDR0001occmax,     PREoccmax, POSToccmax, DIFFoccmax,     Pgamma_max,      maxPow, maxF, maxWF, maxWPow
%%
filename = strcat(savemegto, '/', subj, '_Lcmv_source_spectra.mat');
save (filename, 'readme', 'EPOCH', 'PREoccsel', 'POSToccsel', 'DIFFoccsel', 'PREoccmax', 'POSToccmax', 'DIFFoccmax', 'Nvoxels', 'FDR05occsel', 'FDR0001occsel', 'FDR05occmax', 'FDR0001occmax',...
          'Freq', 'Pgamma_max', 'Pgamma_occmax', 'Pgammarange_sel',...
          'selPow', 'selF', 'selWF', 'selWPow', 'maxPow', 'maxF', 'maxWF', 'maxWPow', 'MAX_ROI', 'MAX_coord_occ', 'MAX_coord','gammamax', 'SORTED', 'SORTEDocc',...
          'binsfreqssel', 'binsfreqsmax', 'Spectrum_pre', 'Spectrum_post', 'MAX_VOXEL_OCC_NUMBER', 'MAX_VOXEL_NUMBER', 'NVOXELS_ID', ...
          'selWMAXF', 'selWMAXPow', 'Pgammarange_sel_MAXselection');


%%
exportToPPTX('saveandclose', PPTXname);




