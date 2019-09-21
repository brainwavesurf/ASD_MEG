%%
clear;
close all;
clc;
screensize = get( groot, 'Screensize' );

tapsmofrq  = 2;
megfolder = strcat( 'meg_sensors_tapsmofrq', num2str(tapsmofrq));
alpharange = [7 14];

% NB: add to path: 
fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(path, fieldtripfolder)
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);
mainpath = '/home/a_shishkina/data/KI/';
realdatapath = '/home/a_shishkina/data/KI/SUBJECTS/';
savepath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';
%%
%add list of subjects:
%SUBJ = ['0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0136'; '0137'; '0138'; '0139'; '0140'; '0141'; '0158'; '0159'; '0160'; '0161'; '0162'; '0163'; '0164'; '0178'; '0179'; '0253'; '0254'; '0255'; '0256'; '0257'; '0259'; '0273'; '0274'; '0275'; '0276'; '0277'; '0346'; '0347'; '0348'; '0350'; '0351'; '0357'; '0358'; '0378'; '0380'; '0381'; '0382'; '0383'; '0384'; '0385']; 
SUBJ = ['0076'];
%%
for s=1: size (SUBJ,1)
    close all
    subj = SUBJ (s,:); 
    savemegto = strcat(savepath, subj);
    epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
    
    %% Load unfiltered epochs. If you need to filter the data (e.g. for LCMV), import raw, not epochs (commented code below).
    ep_fiff_file = strcat(epofolder, subj, '-noerror-lagcorrected-epo.fif');
    hdr = ft_read_header(ep_fiff_file);
    
    
    cfg = [];  
    cfg.dataset = ep_fiff_file;
    cfg.channel={'meg'};
    epochs = ft_preprocessing(cfg);
    
    %%  group preceding events in order to select the epochs later on
    load ([ savemegto, '/', subj, '_info.mat'])
    ev1 = find(allinfo.prev_stim_type==2); % i.e. epochs following Slow (2)
    %ev2 = find(allinfo.prev_stim_type==4); %... medium
    ev2 = find(allinfo.prev_stim_type==8); % ... fast
    prev = {ev1, ev2};
   
    
    %% posterior sensors
    %ch_post = {'MEG1932',  'MEG1922', 'MEG2042',  'MEG2032',  'MEG2112', 'MEG2122',  'MEG2342', 'MEG2332',  'MEG1732', 'MEG1942', 'MEG1912', 'MEG2012', 'MEG2022', 'MEG2312', 'MEG2322', 'MEG2512',...
    %      'MEG1933',  'MEG1923', 'MEG2043',  'MEG2033',  'MEG2113', 'MEG2123',  'MEG2343', 'MEG2333',  'MEG1733', 'MEG1943', 'MEG1913', 'MEG2013', 'MEG2023', 'MEG2313', 'MEG2323', 'MEG2513'};

    %%    
    for con=1:2 % for conditions
        %select epochs according to preceding trials
        
        cfg = [];
        cfg.trials = prev{con};
        cfg.channel = epochs.label;
        cfg.bpfreq = alpharange;
        dataPre = ft_preprocessing(cfg, epochs);
        dataPre_alpha{con} = dataPre;
        
        cfg = [];
        cfg.latency = [-0.8 0.0];
        cfg.covariance = 'yes';
        avg = ft_timelockanalysis(cfg, dataPre_alpha{con});
        covar{con} = avg.cov;
        
        cfg=[];
        cfg.frequency = alpharange; % Hz
        dataPre_tr_avg = ft_selectdata (cfg, dataPre_alpha{con});
        dataPre_avg{con} = dataPre_tr_avg; 
       
        
    end
    %% do eigenvalue decomposition for afterFast and afterSlow
    
    [eigvec, eigval] = eig(covar{2}, covar{2} + covar{1});
    
    spatial_filt = eigvec; % spatial filter / mixing matrix
    spatial_filt_fast = eigvec(:,1); % eigenvector corresponding to min variance for fast condition 
    spatial_filt_slow = eigvec(:,end); % eigenvector corresponding to min variance for slow condition
    
    spatial_patt = (inv(eigvec))'; % spatial pattern / unmixing matrix
    spatial_patt_fast = spatial_patt(:,1); 
    spatial_patt_slow = spatial_patt(:,end);
    
%     %% Whitening process
%     White = sqrt(inv(diag(eigval))) * eigvec'; %Find Whitening Transformation Matrix - Ramoser Equation 
%     
%     S{3} = White * covar{3} * White'; %Whiten Data Using Whiting Transform - Ramoser Equation 
%     S{1} = White * covar{1} * White';
%     
%     [V,D] = eig(S{3},S{1}); %generalized eigenvectors/values
%     % Simultanous diagonalization. Should be equivalent to [B,D]=eig(S{1});
%     [D,ind]=sort(diag(D));
%     V=V(:,ind)
%     
%     %Resulting Projection Matrix-these are the spatial filter coefficients
%     spatial_filt = V'*White;
%     spatial_patt = (inv(spatial_filt))'; %spatial patterns
%     
    %% CSP_data
    %ch x time matrix
    Xfast = dataPre_avg{2}.trial;
    Xslow = dataPre_avg{1}.trial;
    
    Xfast_CSP = spatial_filt_fast' * Xfast{1};
    Xfast_CSP_trial = {Xfast_CSP};
    Xslow_CSP = spatial_filt_slow' * Xslow{1};
    Xslow_CSP_trial = {Xslow_CSP};
    
    %% load atlas, template grid, common for all subjects
    %  load atlas
    atlas = ft_read_atlas( strcat (fieldtripfolder, '/template/atlas/aal/ROI_MNI_V4.nii') ); 
    atlas = ft_convert_units(atlas,'cm');% assure that atlas and template_grid are expressed in the %same units

    %  template MRI
    templatefile = strcat (fieldtripfolder, '/external/spm8/templates/T1.nii'); 
    template_mri = ft_read_mri(templatefile);
    template_mri.coordsys = 'mni';

    %% 
    % load template grid. Use the same template grid, that was used for the
    % warped grid construction!
    load (strcat(mainpath, 'FT_beamformer/', subj, '/mri_linwarp_6mm_brthr0.5/', subj, '_template_grid.mat'));
    % load leadfield    
    load (strcat(mainpath, 'FT_beamformer/', subj, '/mri_linwarp_6mm_brthr0.5/', subj, '_grid_MNI_lf.mat'));
    %%
    cfg = [];
    cfg.atlas = atlas;
    cfg.roi = atlas.tissuelabel;  % here you can also specify a single label, i.e. single ROI
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
    %% perform source analysis
    sourcemodel = grid_MNI_lf;
    sourcemodel.inside = reshape(template_grid.inside, [1, template_grid.dim(1)*template_grid.dim(2)*template_grid.dim(3)]);
    
    %%  
    % Multiply filters with the data and organize into FieldTrip sensable data structure

    virtsens_fast = []; 
    virtsens_slow = [];
    
    for i=1:length(dataPre_avg{2}.trial)
        virtsens_fast.trial{i} = spatial_filt*dataPre_avg{2}.trial{i}; 
    end
        
    virtsens_fast.time = dataPre_avg{2}.time;
    virtsens_fast.fsample = dataPre_avg{2}.fsample;
    virtsens_fast.label= dataPre_avg{2}.label{1};; %(occid(isnotempt)))'
    virtsens_fast.pos = template_grid.pos(find(sourcemodel.inside),:);   
    virtsens_fast.grad = dataPre_avg{2}.grad;
    %virtsens_pre.sampleinfo = data2.sampleinfo;

    for i=1:length(dataPre_avg{1}.trial)
        virtsens_slow.trial{i} = spatial_filt*dataPre_avg{1}.trial{i}; 
    end
    virtsens_slow.time = dataPre_avg{1}.time;
    virtsens_slow.fsample = dataPre_avg{1}.fsample;
    virtsens_slow.label = dataPre_avg{1}.label{1}; %(occid(isnotempt)))'
    virtsens_slow.pos = template_grid.pos(find(sourcemodel.inside),:);   
    virtsens_slow.grad = dataPre_avg{1}.grad;
    %virtsens_post.sampleinfo = data2.sampleinfo;
    %% do spectral analysis
    cfg = [];
    cfg.method    = 'mtmfft';
    cfg.output    = 'pow';
    cfg.foilim    = [5 alpharange(2)];
    cfg.pad = 'nextpow2';
    cfg.tapsmofrq = tapsmofrq;
    cfg.keeptrials = 'yes';
    freq1          = ft_freqanalysis(cfg, virtsens_pre); % keeptrial
    freq2          = ft_freqanalysis(cfg, virtsens_post);
    %% Frequency analysis
    dataPre_avg{2}.trial = Xfast_CSP_trial;
    dataPre_avg{1}.trial = Xslow_CSP_trial;
    
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.output ='pow'; % 'fourier'
    cfg.taper = 'hanning';
    cfg.keeptrials = 'yes';
    cfg.pad = 10; % to [vertually] increase spectral resolution
    cfg.tapsmofrq = tapsmofrq;
    cfg.foilim = alpharange; %freq band of interest
    freqPre_fast = ft_freqanalysis(cfg, virtsens_fast);
    freqPre_slow = ft_freqanalysis(cfg, virtsens_slow);
    
    
    freqPre_fast.powspctrm = freqPre_fast.powspctrm';
    freqPre_slow.powspctrm = freqPre_slow.powspctrm';
        
    MAXmag_fast = max(freqPre_fast.powspctrm(1:3:306));
    MAXmag_slow = max(freqPre_slow.powspctrm(1:3:306));
    MAXgrad_fast = max(freqPre_fast.powspctrm(3:3:204));% for plot, to use the same scale for all events
    MAXgrad_slow = max(freqPre_slow.powspctrm(3:3:204));
    
    %% calculate absolute power differences  between conditions: afterFast - afterSlow
    diff_v3v1=freqPre_fast; % use freqPre{1} structure, but change power for the 'difference'
    diff_v3v1.powspctrm = freqPre_fast.powspctrm  - freqPre_slow.powspctrm;
    %%
    pos_fig1 = [10 10 screensize(3)/10*9 screensize(4)/3];    
    hh=figure('Position',pos_fig1);

    cfg=[];
    cfg.layout = 'neuromag306planar.lay'; % plot only gradeometers
    cfg.zlim = [-1*max(MAXgrad_fast), max(MAXgrad_fast)];
    cfg.colorbar = 'yes';
    subplot(1,4,1)
    ft_topoplotER(cfg,freqPre_slow); title('following Slow');
    subplot(1,4,2)
    ft_topoplotER(cfg,freqPre_fast); title('following Fast');
    subplot(1,4,3); ft_topoplotER(cfg,diff_v3v1); title('alpha power Fast-Slow');
    
    suptitle ('Gradiometers Alpha power [-0.8 to 0 s], according to the type of the preceeding event')
    
end

