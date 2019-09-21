% ft_step2_VirtualCh_lcmv
% cfg.method='lcmv' beamformer.subj
%
% Calculate spectral power in virtual sensors
% Calculate changes in the gamma range 

%For good solution this program uses the following:
% 1. Coregistration should be perfect: if source is in cerebellum, check the
% coregistration!
% 2. Use common covariance for all data (exclude 'active ERP window' from
% analysis [-0.9 -0.1]; [0.4 1.2])
% 3. Use DFT filter on continious signal before epoching to remove 50 Hz
% artefact

%%
clear
close all

%SUBJ = [ 'A001'; 'A002';  'A003'; 'A004'; 'A005'; 'A006'; 'A007'; 'A008'; 'A009'; 'A010'; 'A011'; 'A012'; 'A013'; 'A014'; 'A015'; 'A016'; 'A017']; %'A003';'A001';
%SUBJ = [  'A003'; 'A004'; 'A005'; 'A006'; 'A007'; 'A008'; 'A009'; 'A010']; 
%SUBJ = [ ];
%SUBJ = [ ];



screensize = get( groot, 'Screensize' );

fieldtripfolder = '/mnt/home/a_shishkina/fieldtrip/';
%fieldtripfolder = '/Applications/fieldtrip-20161020/';

path(path, fieldtripfolder)
ft_defaults;

path( strcat (fieldtripfolder, 'external/mne/'), path);
path(path,'/mnt/home/a_shishkina/externals/pptx/');
path(path,'/mnt/home/a_shishkina/projects/asd_meg/0101/MEGdata/FT/'); %for fdr_bh.m

%realdatapath = '/Users/mtw/MEG/Karolinska/MEG_SUBJECTS/'; % 
realdatapath = '/mnt/home/a_shishkina/data/KI/SUBJECTS/'; 
%DataPath = '/Users/mtw/MEG/Karolinska/FT_beamformer/';
DataPath = '/mnt/home/a_shishkina/data/KI/FT_beamformer/';

tapsmofrq = 2; %25
gridres = 6; % 6 mm grid
lambda='5%';
%tapsmofrqnew = 5;

alpharange = [5 20];
alphamax = [7 14];

Nvoxels = 25; % number of voxels [closest to the max] to average.

filt = 'yes';
lpfreq=alpharange(2)+5;
hpfreq=[];

% colors for spectra


color{1}{1}=':b'; color{1}{2}='--b'; color{1}{3}='b';
color{2}{1}=':g'; color{2}{2}='--g'; color{2}{3}='g';
color{3}{1}=':r'; color{3}{2}='--r'; color{3}{3}='r';


megfolder = strcat( '/meg', num2str(gridres), 'mm_linwarp_tapsmofrq', num2str(tapsmofrq), '_lcmv_alpha_', '/');

SUBJ = ['0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0136'; '0137'; '0138'; '0139'; '0140'; '0141'; '0158'; '0159'; '0160'; '0161'; '0162'; '0163'; '0164'; '0178'; '0179'; '0253'; '0254'; '0255'; '0256'; '0257'; '0259'; '0273'; '0274'; '0275'; '0276'; '0277'; '0346'; '0347'; '0348'; '0350'; '0351'; '0357'; '0358'; '0378'; '0380'; '0381'; '0382'; '0383'; '0384'; '0385']; % 
                               

% Segmentations with different brain thresholds:
list05 = ['0380';  '0382';  '0383';  '0384';  '0385'; '0076';  '0102'; '0107'; '0139'; '0140'; '0141'; '0160'; '0161'; '0163'; '0178'; '0254'; '0259'; '0273'; '0274'; '0275';  '0277'; '0346'; '0347'; '0348'; '0350'; '0358' ]; %
list09 = [ '0103'; '0104'; '0106'; '0136'; '0137'; '0138'; '0158'; '0159';'0164'; '0179'; '0255'; '0257'; '0351'; '0357'; '0378'];
list098 = [ '0105'; '0162'; '0253'; '0256';'0381'] ;
list03 = ['0276'; '0101'];
list01 = ['0276']; 
%list05 = [ 'A001';'A002'; 'A003'; 'A004'; 'A005'; 'A006'; 'A007'; 'A008'; 'A009'; 'A010'; 'A011'; 'A012'; 'A013'; 'A014'; 'A015'; 'A016'; 'A017'];

%% load atlas, template grid, common for all subjects

subj = '0101';
%  load atlas
atlas = ft_read_atlas( strcat (fieldtripfolder, '/template/atlas/aal/ROI_MNI_V4.nii') ); 
atlas = ft_convert_units(atlas,'cm');% assure that atlas and template_grid are expressed in the %same units

%  template MRI
templatefile = strcat (fieldtripfolder, '/external/spm8/templates/T1.nii'); 
template_mri = ft_read_mri(templatefile);
template_mri.coordsys = 'mni';

% % 
% load template grid. Use the same template grid, that was used for the
% warped grid construction!
load ('/home/a_shishkina/data/KI/FT_beamformer/0101/mri_linwarp_6mm_brthr0.3/0101_template_grid.mat');

% load ( strcat(DataPath , subj, mrifolder, subj, '_template_grid.mat') );
%

cfg = [];
cfg.atlas      = atlas;
cfg.roi        = atlas.tissuelabel;  % here you can also specify a single label, i.e. single ROI
cfg.inputcoord = 'mni';
mask          = ft_volumelookup(cfg, template_grid);

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
occid2   = find(strncmpi(ROI,'Calcarine',9));  %  indice v1
occid3   = find(strncmpi(ROI,'Cuneus',6));  %  indice
occid4   = find(strncmpi(ROI,'Lingual',7));  %  indice
occid5   = find(strncmpi(ROI,'Precuneus', 9));  
occid    = sort([occid1, occid2, occid3, occid4, occid5]);        
OCC      = ROI(occid);  % label

mask_occ = zeros(size(template_grid.pos,1), 1); mask_occ(ids(occid)) = 1;

%% across subjects
%s=1: size (SUBJ,1)

for s=1: size (1,1)
    close all
    subj = SUBJ (s,:);

    mkdir (strcat(DataPath, subj, megfolder));
    savemegto = strcat(DataPath, subj, megfolder);
    % Save spectra and diff to
    filename = strcat(savemegto, subj, '_', '_Lcmv_source_spectra.mat');

    %mrifolder = strcat( '/mri', num2str(gridres), 'mm', '/');

    if strmatch(subj,list05) 
        mrifolder =strcat('/mri_linwarp_', num2str(gridres), 'mm_', 'brthr',  '0.5/');
    elseif strmatch(subj,list09) 
        mrifolder =strcat( '/mri_linwarp_', num2str(gridres), 'mm_', 'brthr',  '0.9/');
    elseif strmatch(subj,list098) 
        mrifolder =strcat( '/mri_linwarp_', num2str(gridres), 'mm_', 'brthr',  '0.98/');
    elseif strmatch(subj,list03) 
        mrifolder =strcat('/mri_linwarp_', num2str(gridres), 'mm_', 'brthr',  '0.3/');
    elseif strmatch(subj,list01) 
        mrifolder =strcat('/mri_linwarp_', num2str(gridres), 'mm_', 'brthr',  '0.1/');
    else
    %         print ('Subject is not in the list!');
    %         return
        mrifolder =strcat('/mri_linwarp_', num2str(gridres), 'mm_', 'brthr',  '0.5/');

    end

    %% start PPTX report
    exportToPPTX('new');
    PPTXname  = strcat(savemegto, subj,'_', '_lcmv_report');

    %% subject specific Loads
    % load leadfield / % source model    
    load ( strcat(DataPath , subj, mrifolder, subj, '_grid_MNI_lf.mat') );
    % head model: individ_hdm_vol
    load ( strcat(DataPath , subj, mrifolder, subj, '_individ_hdm_vol.mat' )); 
    % DICS results
    %     load (strcat(DataPath , subj, megfolder, subj, '_DISCbroad_source_maxT.mat'));
    % load individual mri
    load ( strcat(DataPath , subj, mrifolder, subj, '_mri_orig_realigned.mat' )); 

    %% read clean and delay corrected events exported from mne .py (exported with events_raw_resample_adults_4vel.py)
    load (strcat (realdatapath, subj, '/ICA_nonotch_crop/', subj, '_clean_events.mat'));

    %% Load data epochs 
    %     ep_fiff_file = strcat(realdatapath, subj, '/epochs/', subj, '-noerror-lagcorrected-epo.fif');
    %     hdr = ft_read_header(ep_fiff_file);

    %     %% find good epochs
    %     C = strsplit(hdr.orig.epochs.drop_log,'], ');
    %     find1=strcmp(C, '[[');
    %     find2=strcmp(C, '[');
    %     find3=strcmp(C, '[]]');
    %     ind = find(find1+find2+find3);
    %     % read events (saved by '/Users/mtw/MEG/Scripts/Karolinska/PY/Sensors/resample_raw.py' )
    %     events = load (strcat (realdatapath, subj, '/ICA_nonotch_crop/', subj, '_events.mat'));
    %     events = events.events( ind, :);
    %% Load raw data
    fiff_file = strcat(realdatapath, subj, '/ICA_nonotch_crop/', subj, '_rings_ICA_raw.fif');
    hdrraw = ft_read_header(fiff_file);
    pre = -1.0* hdrraw.Fs ;
    post = 1.2* hdrraw.Fs ;
    trl=[];
    first= round(cast(hdrraw.orig.raw.first_samp, 'double'));
    events(:,1) =  events(:,1)-first;
    
    for i=1:size (events,1)
         trl(i, 1)=(events(i,1)+pre) ; 
         trl(i, 2)=(events(i,1)+post) ; 
         trl(i, 3)= -1.0*hdrraw.Fs ; % offset
         trl(i, 4) = events(i,3); % stimulus_value;
    end

    %% extract data and epochs from the raw
    cfg = [];
    cfg.trl=trl;
    cfg.channel     = 'meg';
    cfg.dftfilter   = 'yes';
    cfg.dftfreq     = [50 100];
    cfg.demean = 'yes';
    if strcmp(filt,'yes')
       cfg.lpfilter  = 'yes';
       cfg.lpfreq    = lpfreq;
       %cfg.hpfilter  = 'yes';
       %cfg.hpfreq    = hpfreq;
    end
    cfg.dataset = fiff_file;
    cfg    = ft_definetrial(cfg);
    epochs = ft_preprocessing(cfg);

    %% Plot average
    cfg = [];
    cfg.channel=epochs.label;
    avg = ft_timelockanalysis(cfg,epochs);

    hh=figure;
    ttt = find(avg.time>-0.1 & avg.time<0.5);
    plot (avg.time(ttt), avg.avg(:, ttt));
    title ('Sensor averages');

    exportToPPTX('addslide'); % slide with distributions
    exportToPPTX('addpicture', hh, 'Position', [0.5,0.5,9,6]);


    %%  select epochs according to events

    ev1 = find(events(:,3)==2);
    ev2 = find(events(:,3)==4);
    ev3 = find(events(:,3)==8);
    EV = {ev1, ev2, ev3};

    %%
    cfg = [];
    cfg.channel           = epochs.label;
    cfg.covariance        = 'yes';
    %cfg.covariancewindow  =  [-0.9 -0.1; 0.3 1.1];%'all'; %[-0.9 -0.1]; %
    cfg.covariancewindow  =  'all'; %[-0.9  1.2];%'all'; %[-0.9 -0.1]; %
    cfg.vartrllength      = 2;
    total_avg = ft_timelockanalysis(cfg, epochs);

    %covariance.wind = cfg.covariancewindow;
    %covariance.epochs = 'all conditions';

    %% perform source analysis
    % sourcemodel = load ( strcat(DataPath , subj, '/mri/', subj, '_subj_warped_grid.mat') );
    % sourcemodel = sourcemodel.subj_warped_grid;
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
    cfg.normalize = 'yes';
    source_total=ft_sourceanalysis(cfg, total_avg);
    source_total.pos = template_grid.pos;


    %% 
    
    for con = 1:3 % for conditions 
    %% select epochs
        cfg = [];
        cfg.trials =  EV{con};        
        EPO = ft_selectdata(cfg, epochs);

        cfg = [];
        cfg.channel=EPO.label;
        avg = ft_timelockanalysis(cfg,EPO);
        AVG{con}=avg;

        % subtract
        %for ttt=1:size (EPO.trial,2)
        %    EPO.trial{ttt}=EPO.trial{ttt}-avg.avg;
        %end

        %% select pre and poststim   
        cfg=[];
        cfg.latency     = [-0.9 0];
        data2_pre = ft_selectdata(cfg, EPO);

        cfg.latency     = [0.3 1.2];
        data2_post = ft_selectdata(cfg, EPO);

        %% This is for cfg.lcmv.fixedori='yes'; 
        % Multiply filters with the data and organize into FieldTrip sensable data structure
        spatialfilter=cat(1,source_total.avg.filter{:});  % all voxels inside
        %spatialfilter_post=cat(1,source_total.avg.filter{:});  % all voxels inside

        virtsens_pre=[]; virtsens_post=[];
        for i=1:length(data2_pre.trial)
            virtsens_pre.trial{i}=spatialfilter*data2_pre.trial{i};
            virtsens_post.trial{i}=spatialfilter*data2_post.trial{i};
        end
        virtsens_pre.time=data2_pre.time;
        virtsens_pre.fsample=data2_pre.fsample;
        virtsens_pre.label= cellstr(string(find(sourcemodel.inside))); %(occid(isnotempt)))'
        virtsens_pre.pos = template_grid.pos(find(sourcemodel.inside),:);   
        virtsens_pre.grad =data2_pre.grad;
        %virtsens_pre.sampleinfo = data2.sampleinfo;

        virtsens_post.time=data2_post.time;
        virtsens_post.fsample=data2_post.fsample;
        virtsens_post.label= cellstr(string(find(sourcemodel.inside))); %(occid(isnotempt)))'
        virtsens_post.pos = template_grid.pos(find(sourcemodel.inside),:);   
        virtsens_post.grad =data2_post.grad;
        %virtsens_post.sampleinfo = data2.sampleinfo;

        %%
        % % aaa=cat(3,virtsens.trial{:});
        % % ave_source=squeeze(mean(aaa,3));
        % % % % aa = std(reshape (aaa, [size(aaa,1), size(aaa,3)*size(aaa,2)]),0,2); 
        % % hh=figure;
        % % subplot(2,1,1);
        % % plot (avg.time(ttt), avg.avg(:, ttt));
        % % title ('sensor averages')
        % % 
        % % subplot(2,1,2);
        % % plot (avg.time(ttt), ave_source(:, ttt));
        % % title ('Source averages')
        % % exportToPPTX('addslide'); % slide with distributions
        % % exportToPPTX('addpicture', hh, 'Position', [0.5,0.5,9,6]);

        %% This is for cfg.lcmv.fixedori='no'; 
        % % % %% Multiply filters with the data and organize into FieldTrip sensable data structure
        % % % spatialfilter=cat(3,source.avg.filter{:});  % all voxels inside
        % % % virtsens=[];
        % % % for i=1:length(data2_pre.trial)
        % % %     i
        % % %     for j =1:size(spatialfilter,3)
        % % %         tmp.trial{i,j}=spatialfilter(:,:,j)*data2_pre.trial{i};  % virtual ch {trial, voxel}
        % % %     end
        % % % end
        % % % %% construct a single virtual channel in the maximum power orientation
        % % % 
        % % % virtsens=[];
        % % % for j=1:size(spatialfilter,3)
        % % %     timeseries = cat(2,  tmp.trial{:,j});
        % % %     [u, s, v] = svd(timeseries, 'econ');
        % % %     j
        % % %     for i=1:length(data2.trial)
        % % %       virtsens.trial{i}(j,:) = u(:,1)' * squeeze(spatialfilter(:,:,j)) * data2_pre.trial{i};  % filter dim is  3 x 306
        % % %     end
        % % % 
        % % % end
        % % % %%
        % % % virtsens.trial = virtuasens.trial;
        % % % virtsens.time=data2.time;
        % % % virtsens.fsample=data2.fsample;
        % % % virtsens.label= cellstr(string(find(sourcemodel.inside))); %(occid(isnotempt)))'
        % % % virtsens.pos = sourcemodel.pos(find(sourcemodel.inside),:);   
        % % % virtsens.grad =data2.grad;
        % % % virtsens.sampleinfo = data2.sampleinfo;

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
        cfg.keeptrials = 'no';
        fd_pre         = ft_freqdescriptives(cfg, freq1); % average
        fd_post        = ft_freqdescriptives(cfg, freq2);

        % %     %% Spectrum in the maximally modulated voxel plus surround
        % %     % load saved: maxindOCC from 
        % %     loadfilename = strcat(DataPath, subj, '/',  'meg6mm_linwarp_tapsmofrq5_individWF_maxchange_static/', subj, '_DISCbroad_source_maxT.mat');
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
        f1=alphamax(1); f2=alphamax(2);
        alpha_ind = find(fd_pre.freq>=floor(f1) & fd_pre.freq<=ceil(f2) ); % alpha frequencues
        alpha_wholerange_ind = find(fd_pre.freq>=floor(alpharange(1)) & fd_pre.freq<=ceil(alpharange(2)) ); % alpha frequencues
        inside_id = find(template_grid.inside);
        %DIFFalpha=   mean( fd_post.powspctrm(:, alpha_ind)-fd_pre.powspctrm(:, alpha_ind)./fd_pre.powspctrm(:, alpha_ind) ,2); % average diff in alpha band
        DIFFalpha=   (mean( fd_post.powspctrm(:, alpha_ind),2)-mean(fd_pre.powspctrm(:, alpha_ind),2))./mean(fd_pre.powspctrm(:, alpha_ind) ,2); % average diff in alpha band

        Spectrum_pre{con} = fd_pre;
        Spectrum_post{con} = fd_post;


        % plot DIFFalpha
        SourseDIFF = source_total;
        SourseDIFF.pow = DIFFalpha;
        SourseDIFF = rmfield(SourseDIFF,'avg')
        SourseDIFF.pow = NaN(size(SourseDIFF.inside))';
        SourseDIFF.pow(find(SourseDIFF.inside)) = DIFFalpha';
        SourseDIFF.coordsys = 'spm'
        
        savefolder = '/mnt/home/a_shishkina/data/KI/Save/'  
        save(savefolder, 'SourseDIFF')
        
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
        cfg.location = 'min';
        %cfg.opacitylim    = 'zeromax'; 
        cfg.opacitymap    = 'rampup'; 
        cfg.atlas=atlas;
        %cfg.roi = {'Occipital_Sup_R', 'Occipital_Sup_L',  'Occipital_Inf_R',  'Occipital_Inf_L', 'Cuneus_R', 'Cuneus_L', 'Precuneus_R' , 'Precuneus_L', 'Occipital_Mid_R' , 'Occipital_Mid_L', 'Calcarine_R', 'Calcarine_L', 'Lingual_R' , 'Lingual_L'      };
        ft_sourceplot(cfg, SourseDIFF_Int );

        % save distribution of the ratio:
        exportToPPTX('addslide'); % slide with distributions
        exportToPPTX('addpicture', gcf, 'Position', [0.3,1,4,3]);
        exportToPPTX('addtext', strcat('lcmv beamformer; Subject:', subj, ',total and occipital,  condition V', num2str(con), ': (post-pre)/pre in the alpha range: ', num2str(alphamax(1)), '-', num2str(alphamax(2)),'Hz, filter=', filt  ), 'FontSize', 14, 'Position', [0.5, 0.5, 10,1]);

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
        cfg.location = 'min';
        %cfg.opacitylim    = 'zeromax'; 
        cfg.opacitymap    = 'rampup'; 
        cfg.atlas=atlas;
        ft_sourceplot(cfg, SourseDIFF_Int );
        exportToPPTX('addpicture', gcf, 'Position', [4.5,1,4,3]);


        %% Save positions and probability for maximally modulated voxel, whole brain
        [DIFFmax,x] = max(DIFFalpha);    

        A = [DIFFalpha, inside_id, [1:length(inside_id)]', virtsens_post.pos ];
        B=sortrows(A,1); 
        max_voxel_inside = B(1,3); %index inside
        max_voxel_pos = virtsens_post.pos(max_voxel_inside,:); % Max alpha SUPPRESSION!
        MAX_VOXEL_NUMBER{con} = max_voxel_inside;

        %%%% for plot 
        n=1:Nvoxels;
        max_voxels_inside = B(n,3);
        max_voxels_pos = virtsens_post.pos(max_voxels_inside,:);
        MAX_ROI{con}= ROI(B(n,3));
        MAX_coord{con}= source_total.pos(B(n,2),:);
        %%%%%


        % probability of in the max modulated voxel
        pre =  squeeze(mean(freq1.powspctrm(:, max_voxel_inside, alpha_wholerange_ind),2)); % tr, voxel_inside, freq; / average over max_voxels_inside and alpha_ind
        post = squeeze(mean(freq2.powspctrm(:, max_voxel_inside, alpha_wholerange_ind),2));
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

        % p for averaged alphamax range
        pre =  squeeze(mean(mean(freq1.powspctrm(:, max_voxel_inside, alpha_ind),3),2)); % tr, voxel_inside, freq; / average over max_voxels_inside and alpha_ind
        post = squeeze(mean(mean(freq2.powspctrm(:, max_voxel_inside, alpha_ind),3),2));
        [Palpha_max{con}, h, statis] = signrank(pre,post); % P in the alpha range in the max voxel


        %%  Save positions and probability for maximally modulated voxel, occipital region
        % ids, occid
        Aocc = [DIFFalpha(occid), inside_id(occid), occid', virtsens_post.pos(occid,:) ];
        Bocc =  sortrows(Aocc,1);
        max_voxel_inside_occ = Bocc(1,3);
        max_voxel_pos_occ = virtsens_post.pos(max_voxel_inside_occ,:); % Max alpha SUPPRESSION!
        MAX_VOXEL_OCC_NUMBER{con} = max_voxel_inside_occ;

        %%%% for plot 
        n=1:Nvoxels;
        max_voxels_inside_occ = Bocc(n,3);
        max_voxels_pos_occ = virtsens_post.pos(max_voxels_inside_occ,:);
        MAX_ROIocc{con}= ROI(Bocc(n,3));
        MAX_coord_occ{con}= source_total.pos(Bocc(n,2),:);
        %%%%%

        % probability in the max modulated occ voxels
        pre_occ =  squeeze(freq1.powspctrm(:, max_voxel_inside_occ, alpha_wholerange_ind) ); % tr, voxel_inside, freq; / average over max_voxels_inside and alpha_ind
        post_occ = squeeze(freq2.powspctrm(:, max_voxel_inside_occ, alpha_wholerange_ind) );
        for j=1: size(pre_occ,2) % for all freq
           [h, p] = ttest(squeeze(pre_occ(:,j)),squeeze(post_occ(:,j)));
           PValues(j)=p;
        end
        [FDR05occmax{con}.h, FDR05occmax{con}.crit_p, FDR05occmax{con}.adj_p]=fdr_bh(PValues ,0.0001 ,'pdep','yes');
        [FDR0001occmax{con}.h, FDR0001occmax{con}.crit_p, FDR0001occmax{con}.adj_p]=fdr_bh(PValues ,0.0001 ,'pdep','yes');
        maxvoxel_pre_occ{con} = fd_pre.powspctrm(max_voxel_inside_occ,:);
        maxvoxel_post_occ{con} = fd_post.powspctrm(max_voxel_inside_occ,:);
        MAX_ROIocc{con}= ROI(Bocc(1,3));
        MAX_coord_occ{con}= source_total.pos(Bocc(1,2),:);

       % p for averaged alphamax range
        pre =  squeeze(mean(freq1.powspctrm(:, max_voxel_inside_occ, alpha_ind),3)); % tr, voxel_inside, freq; / average over max_voxels_inside and alpha_ind
        post = squeeze(mean(freq2.powspctrm(:, max_voxel_inside_occ, alpha_ind),3));
        [Palpha_occmax{con}, h, statis] = signrank(pre,post); % P in the alphamax range in the occipital max voxel, for averaged alphamax range


        PREoccmax{con}  = squeeze(fd_pre.powspctrm(max_voxel_inside_occ, :));
        POSToccmax{con} = squeeze(fd_post.powspctrm(max_voxel_inside_occ, :));
        DIFFoccmax{con} = (POSToccmax{con}-PREoccmax{con})./PREoccmax{con};

        % wF and power for DIFFoccmax
        alphaind = find(fd_post.freq>=alphamax(1) &  fd_post.freq<=alphamax(2));
        [maxPow{con},  ind] = min(DIFFoccmax{con}(alphaind)); % SUPPRESSION!
        maxF{con} = fd_post.freq(alphaind(ind));
        inds = find(DIFFoccmax{con}(alphaind) <=2/3*maxPow{con});
        % if frequencies are too far from the maximum SUPPRESSION (more than 10 Hz) we exclude them!
        excl = find( (fd_post.freq(alphaind(inds))-maxF{con}) >10);
        inds(excl) =[];
        binsfreqsmax{con} = fd_post.freq(alphaind(inds));

        maxWF{con} = sum(fd_post.freq(alphaind(inds)).*DIFFoccmax{con}(alphaind(inds) ))/sum(DIFFoccmax{con}(alphaind(inds) ));
        maxWPow{con} = mean(DIFFoccmax{con}(alphaind(inds)) );


        %FDR05occmax, FDR0001occmax,    PREoccmax, POSToccmax, DIFFoccmax,     Palpha_max,      maxPow, maxF, maxWF, maxWPow

        %% average power around maxmodulated OCC voxel: Nvoxels averaged
        % calculate distance from the max voxel  MAX_voxel_inside_occ  DIFFalpha
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
        pre_occ =  squeeze(mean(freq1.powspctrm(:,  Nvoxels_id, alpha_wholerange_ind),2)); % tr, voxel_inside, freq; / average over max_voxels_inside and alpha_ind
        post_occ = squeeze(mean(freq2.powspctrm(:,  Nvoxels_id, alpha_wholerange_ind),2));
        for j=1: size(pre_occ,2) % for all freq
           [h, p] = ttest(squeeze(pre_occ(:,j)),squeeze(post_occ(:,j)));
           PValues(j)=p;
        end
        [FDR05occsel{con}.h, FDR05occsel{con}.crit_p, FDR05occsel{con}.adj_p]=fdr_bh(PValues ,0.0001 ,'pdep','yes');
        [FDR0001occsel{con}.h, FDR0001occsel{con}.crit_p, FDR0001occsel{con}.adj_p]=fdr_bh(PValues ,0.0001 ,'pdep','yes');
        selvoxel_pre_occ{con}  = mean(fd_pre.powspctrm(Nvoxels_id,:),1);
        selvoxel_post_occ{con} = mean(fd_post.powspctrm(Nvoxels_id,:),1);

        % p for averaged alphamax range
        pre =  squeeze(mean(mean(freq1.powspctrm(:, Nvoxels_id, alpha_ind),3),2)); % tr, voxel_inside, freq; / average over max_voxels_inside and alpha_ind
        post = squeeze(mean(mean(freq2.powspctrm(:, Nvoxels_id, alpha_ind),3),2));
        [Palpharange_sel{con}, h, statis] = signrank(pre,post); % P in the alpha range in the occipital Nvoxels

        PREoccsel{con}  = mean(fd_pre.powspctrm(Nvoxels_id, :),1);
        POSToccsel{con} = mean(fd_post.powspctrm(Nvoxels_id, :),1);
        DIFFoccsel{con} = (POSToccsel{con}-PREoccsel{con})./PREoccsel{con};

        % wF and power for DIFFoccsel
        [selPow{con},  ind] = min(DIFFoccsel{con}(alphaind));
        selF{con} = fd_post.freq(alphaind(ind));
        inds = find(DIFFoccsel{con}(alphaind) <=2/3*selPow{con});
        % if frequencies are too far from the maximum (more than 20 Hz) we exclude them!
        excl = find( (fd_post.freq(alphaind(inds))-selF{con}) >10);
        inds(excl) =[];
        binsfreqssel{con} = fd_post.freq(alphaind(inds));
        INDS{con}=inds; % alpha band indexes to remember

        selWF{con} = sum(fd_post.freq(alphaind(inds)).*DIFFoccsel{con}(alphaind(inds) ))/sum(DIFFoccsel{con}(alphaind(inds) ));
        selWPow{con} = mean(DIFFoccsel{con}(alphaind(inds)) );


       % FDR05occsel, FDR0001occsel,    PREoccsel, POSToccsel, DIFFoccsel,     Palpha_sel,     selPow, selF, selWF, selWPow


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
            p1=plot (fd_pre.freq(alpha_wholerange_ind), PRE(alpha_wholerange_ind), 'b');  hold on
            p2=plot (fd_pre.freq(alpha_wholerange_ind), POST(alpha_wholerange_ind), 'r'); 
            xlabel('Frequency, Hz')
            ylabel('Power')
            yyaxis right
            ylabel('(pre-post)/pre')
            p3=plot(fd_pre.freq(alpha_wholerange_ind), diff(alpha_wholerange_ind), '--k');
            hold off
            title (  strcat ('Max',  num2str(i), ', ', ROI(B(i,3)), ', coord:', num2str(B(i,4:6), '%4.1f ') ) ) ;        
        end

        legend([p1 p2 p3], {'pre pow', 'post pow', '(pre-post)/pre'}, 'FontSize', 16);
        suptitle(strcat('Subject:', subj, ', pre- and post-stimulus power and (pre-post)/pre ratio in alpha range in condition..', num2str(con), ', 10 maximal voxels in the whole brain') )
        exportToPPTX('addslide'); % slide with distributions
        exportToPPTX('addpicture', hh, 'Position', [0.5,1,9,6]);
        exportToPPTX('addtext', strcat('WHOLE BRAIN:  Subject:', subj, ', Condition:', num2str(con),  ', Filter:', filt, ', alpha range:', num2str(alphamax(1)), '-',...
                     num2str(alphamax(2)), 'Hz ', ', cfg.lcmv.lambda=', lambda, ', tapsmofrq=', num2str(tapsmofrq),  ', Ranksun p for the mean power in 45-110Hz:', num2str(Palpha_max{con})   ),  'FontSize', 14, 'Position', [0.1, 0.1, 10,1]);

      %%  Plot spectra for 10 the maximally modutated voxels, occipital selection
        pos_fig1 = [10 10 screensize(3)/10*9 screensize(4)/10*9];    
        hh=figure('Position',pos_fig1);
        for i=1:N
            i
            POST=mean(fd_post.powspctrm(max_voxels_inside_occ(i),:),1);
            PRE=mean(fd_pre.powspctrm(max_voxels_inside_occ(i),:),1);
            diff = (POST-PRE)./PRE;
            subplot(2,N/2,i);
            p1=plot (fd_pre.freq(alpha_wholerange_ind), PRE(alpha_wholerange_ind), 'b');  hold on
            p2=plot (fd_pre.freq(alpha_wholerange_ind), POST(alpha_wholerange_ind), 'r'); 
            xlabel('Frequency, Hz')
            ylabel('Power')
            yyaxis right
            ylabel('(pre-post)/pre')
            p3=plot(fd_pre.freq(alpha_wholerange_ind), diff(alpha_wholerange_ind), '--k');
            hold off
            title (  strcat ('Max',  num2str(i), ', ', ROI(Bocc(i,3)), ', coord:', num2str(Bocc(i,4:6), '%4.1f ') ) ) ;        
        end
        legend([p1 p2 p3], {'pre', 'post', '(pre-post)/pre'}, 'FontSize', 16);
        suptitle(strcat('Subject:', subj, ', pre- and post-stimulus power and (pre-post)/pre ratio in alpha range in condition..', num2str(con), ', 10 maximal voxels in occipital selection') )

        exportToPPTX('addslide'); % slide with distributions
        exportToPPTX('addpicture', hh, 'Position', [0.5,1,9,6]);
        exportToPPTX('addtext', strcat('OCCIPITAL SELECTION:   Subject:', subj, ', Condition:', num2str(con),  ', Filter:', filt, ', alpha range:', num2str(alphamax(1)), '-',...
        num2str(alphamax(2)), 'Hz ', ', cfg.lcmv.lambda=', lambda, ', tapsmofrq=', num2str(tapsmofrq), ', Minimal FDR-p-level for max-modulated-voxel in alpha range =', num2str(min(FDR05occmax{con}.adj_p)), ', Ranksun p for the mean power in 45-110Hz:', num2str(Palpha_occmax{con})  ), 'FontSize', 14, 'Position', [0.1, 0.1, 10,1]);

        %%
        FREQ1powspctrm{con} = squeeze(mean(freq1.powspctrm(:,:,alpha_ind),3)); %keeptrial
        FREQ2powspctrm{con} = squeeze(mean(freq2.powspctrm(:,:,alpha_ind),3));
    end % end for conditons
    
%     %% power in 25 max voxels based on 100% contrast condition, for contrast 50%
%     if strcmp(experiment, 'contrast') % load the voxel selection from the subject's 100% contrast experiment
%         LC1 = load (strcat('/Volumes/MEGdisk/Moscow_4velocities/Beamformer/', subj, '/meg6mm_linwarp_tapsmofrq5_lcmv_static/', subj, '_static_Lcmv_source_spectra.mat'),  'selWPow'); 
%         LC2 = load (strcat('/Volumes/MEGdisk/Moscow_4velocities/Beamformer/', subj, '/meg6mm_linwarp_tapsmofrq5_lcmv_static/', subj, '_static_Lcmv_source_spectra.mat'), 'NVOXELS_ID');
%         [a,b] = min(cell2mat(LC1.selWPow)); % b is the max condition
%         for con=1:3
%            pre =  squeeze(mean(FREQ1powspctrm{con}(:, NVOXELS_ID{b}),2)); % tr, voxel_inside, freq; / average over max_voxels_inside and alpha_ind
%            post = squeeze(mean(FREQ2powspctrm{con}(:, NVOXELS_ID{b}),2));
%            [Palpharange_sel_MAXselection{con}, h, statis] = signrank(pre,post); % P in the alpha range in the occipital Nvoxels
% 
%            PREoccsel_MAXselection{con}  = mean(Spectrum_pre{con}.powspctrm(NVOXELS_ID{b}, :),1);
%            POSToccsel_MAXselection{con} = mean(Spectrum_post{con}.powspctrm(NVOXELS_ID{b}, :),1);
%            DIFFoccsel_MAXselection{con} = (POSToccsel_MAXselection{con}-PREoccsel_MAXselection{con})./PREoccsel_MAXselection{con};
% 
%            % wFselMAX and power for DIFFoccsel_MAXselection
%            [selMAXPow{con},  ind] = min(DIFFoccsel{con}(alphaind));
%            selMAXF{con} = fd_post.freq(alphaind(ind));
%            inds = find(DIFFoccsel{con}(alphaind) <=2/3*selMAXPow{con});
%            % if frequencies are too far from the maximum (more than 10 Hz) we exclude them!
%            excl = find( abs((fd_post.freq(alphaind(inds))-selMAXF{con})) >10);
% 
%            inds(excl) =[];
%            binsfreqsselMAX{con} = fd_post.freq(alphaind(inds));
% 
%            selWMAXF{con} = sum(fd_post.freq(alphaind(inds)).*DIFFoccsel_MAXselection{con}(alphaind(inds) ))/sum(DIFFoccsel_MAXselection{con}(alphaind(inds) ));
%            selWMAXPow{con} = mean(DIFFoccsel_MAXselection{con}(alphaind(inds)) );

      %% power in 25 max voxels based on the MAX SUPPRESSION condition:  
%     if strcmp(experiment,'static') || strcmp(experiment,'2F2S')
%         [a,b] = min(cell2mat(selWPow)); % b is the max SUPPRESSION condition
%         for con=1:3  
%             pre =  squeeze(mean(FREQ1powspctrm{con}(:, NVOXELS_ID{b}),2)); % tr, voxel_inside, freq; / average over max_voxels_inside and alpha_ind
%             post = squeeze(mean(FREQ2powspctrm{con}(:, NVOXELS_ID{b}),2));
%             [Palpharange_sel_MAXselection{con}, h, statis] = signrank(pre,post); % P in the alpha range in the occipital Nvoxels
% 
%             PREoccsel_MAXselection{con}  = mean(Spectrum_pre{con}.powspctrm(NVOXELS_ID{b}, :),1);
%             POSToccsel_MAXselection{con} = mean(Spectrum_post{con}.powspctrm(NVOXELS_ID{b}, :),1);
%             DIFFoccsel_MAXselection{con} = (POSToccsel_MAXselection{con}-PREoccsel_MAXselection{con})./PREoccsel_MAXselection{con};
% 
%             % wFselMAX and power for DIFFoccsel_MAXselection
%             [selMAXPow{con},  ind] = min(DIFFoccsel{con}(alphaind));
%             selMAXF{con} = fd_post.freq(alphaind(ind));
%             inds = find(DIFFoccsel{con}(alphaind) <=2/3*selMAXPow{con});
%             % if frequencies are too far from the maximum (more than 20 Hz) we exclude them!
%             excl = find( abs((fd_post.freq(alphaind(inds))-selMAXF{con})) >10);
%             inds(excl) =[];
%             binsfreqsselMAX{con} = fd_post.freq(alphaind(inds));
% 
%             selWMAXF{con} = sum(fd_post.freq(alphaind(inds)).*DIFFoccsel_MAXselection{con}(alphaind(inds) ))/sum(DIFFoccsel_MAXselection{con}(alphaind(inds) ));
%             selWMAXPow{con} = mean(DIFFoccsel_MAXselection{con}(alphaind(inds)) );
%          end
%     end
end
%% save figure with spectra in max selection
figure;
plot(fd_pre.freq, PREoccsel{1}, color{1}{1}); hold on; plot(fd_post.freq,POSToccsel{1},color{1}{2});
plot(fd_pre.freq, PREoccsel{2}, color{2}{1}); hold on; plot(fd_post.freq,POSToccsel{2},color{2}{2});
plot(fd_pre.freq, PREoccsel{3}, color{3}{1}); hold on; plot(fd_post.freq,POSToccsel{3},color{3}{2});

xlabel('Frequency, Hz');
ylabel('Power')
yyaxis right
ylabel('(pre-post)/pre')
plot(fd_post.freq, DIFFoccsel{1}, color{1}{3}); 
plot(fd_post.freq, DIFFoccsel{2}, color{2}{3}); 
plot(fd_post.freq, DIFFoccsel{3}, color{3}{3}); 

hold off;
legend( {'SlowPre', 'SlowPost','MediumPre', 'MediumPost', 'FastPre', 'FastPost', 'SlowDiff',  'MediumDiff', 'FastDiff'}, 'FontSize', 14);
title (strcat('LCMV beamformer: selection of..', num2str(Nvoxels), '..voxels around occipital max'), 'FontSize', 14)
exportToPPTX('addslide'); % slide with distributions
exportToPPTX('addpicture', gcf, 'Position', [1,1,7,4.5]);
exportToPPTX('addtext', strcat('Subject:', subj, ', spectral changes in selection of...', num2str(Nvoxels), ' atound max occ. max: lcmv beamformer' ), 'FontSize', 14, 'Position', [0.5, 0.5, 10,1]);
%% save figure with spectra in the max condition selection

% figure;
% plot(fd_pre.freq, PREoccsel_MAXselection{1}, color{1}{1}); hold on; plot(fd_post.freq,POSToccsel_MAXselection{1},color{1}{2});
% plot(fd_pre.freq, PREoccsel_MAXselection{2}, color{2}{1}); hold on; plot(fd_post.freq,POSToccsel_MAXselection{2},color{2}{2});
% plot(fd_pre.freq, PREoccsel_MAXselection{3}, color{3}{1}); hold on; plot(fd_post.freq,POSToccsel_MAXselection{3},color{3}{2});
% plot(fd_pre.freq, PREoccsel_MAXselection{4}, color{4}{1}); hold on; plot(fd_post.freq,POSToccsel_MAXselection{4},color{4}{2});
% 
% xlabel('Frequency, Hz');
% ylabel('Power')
% yyaxis right
% ylabel('(pre-post)/pre')
% plot(fd_post.freq, DIFFoccsel_MAXselection{1}, color{1}{3}); 
% plot(fd_post.freq, DIFFoccsel_MAXselection{2}, color{2}{3}); 
% plot(fd_post.freq, DIFFoccsel_MAXselection{3}, color{3}{3}); 
% plot(fd_post.freq, DIFFoccsel_MAXselection{4}, color{4}{3}); 
% 
% hold off;
% legend( {'StaticPre', 'StaticPost', 'SlowPre', 'SlowPost','MediumPre', 'MediumPost', 'FastPre', 'FastPost', 'StaticDiff', 'SlowDiff',  'MediumDiff', 'FastDiff'}, 'FontSize', 14);
% title (strcat('LCMV beamformer: sel. of..', num2str(Nvoxels), '..voxels around occ. MAXIMAL CONDITION max '), 'FontSize', 14)
% exportToPPTX('addslide'); % slide with distributions
% exportToPPTX('addpicture', gcf, 'Position', [1,1,7,4.5]);
% exportToPPTX('addtext', strcat('Subject:', subj, ', spectral changes in sel. of...', num2str(Nvoxels), ' around occ. MAX CONDITION MAX: lcmv beamformer' ), 'FontSize', 14, 'Position', [0.5, 0.5, 10,1]);


%% save figure with spectra in max 
figure;
plot(fd_pre.freq, PREoccmax{1}, color{1}{1}); hold on; plot(fd_post.freq,POSToccmax{1},color{1}{2});
plot(fd_pre.freq, PREoccmax{2}, color{2}{1}); hold on; plot(fd_post.freq,POSToccmax{2},color{2}{2});
plot(fd_pre.freq, PREoccmax{3}, color{3}{1}); hold on; plot(fd_post.freq,POSToccmax{3},color{3}{2});

xlabel('Frequency, Hz');
ylabel('Power')
yyaxis right
ylabel('(pre-post)/pre')
plot(fd_post.freq, DIFFoccmax{1}, color{1}{3}); 
plot(fd_post.freq, DIFFoccmax{2}, color{2}{3}); 
plot(fd_post.freq, DIFFoccmax{3}, color{3}{3}); 


hold off;
legend( {'SlowPre', 'SlowPost','MediumPre', 'MediumPost', 'FastPre', 'FastPost', 'SlowDiff',  'MediumDiff', 'FastDiff'}, 'FontSize', 14);
title ('LCMV veamformer: max occipital voxel ', 'FontSize', 14)
exportToPPTX('addslide'); % slide with distributions
exportToPPTX('addpicture', gcf, 'Position', [1,1,7,4.5]);
exportToPPTX('addtext', strcat('Subject:', subj, ', spectral changes in  the max occip. voxels: lcmv beamformer' ), 'FontSize', 14, 'Position', [0.5, 0.5, 10,1]);

%% Save spectra and diff
%filename = strcat(savemegto, subj, '_Lcmv_source_spectra.mat');
readme.covariance = 'covatiance window and epochs info'
readme.SORTED = 'sorted by diff=pre-post/pre:  [diff, voxel, voxel_inside, X-position, Y-position, Z-position  ], positions in MNI coordinates'
readme.info = strcat('Subject:', subj, ', Condition:', num2str(con),  ', Filter:', filt, ', alpha range:', num2str(f1), '-',...
             num2str(f2), 'Hz ', ', cfg.lcmv.lambda=', lambda, ', tapsmofrq=', num2str(tapsmofrq)  );
readme.freq = 'frequencies';
readme.FDR = 'false discovery rate corrected p at the [average of] maximally modulated voxels, with 05 and 0001 thresholds';

readme.PREocc  = 'whole power spectrum in the max modulated voxels or selection (of Nxovels), prestimulus';
readme.POSTocc = 'whole power spectrum in the max modulated voxels or selection (of Nxovels), poststimulus';
readme.DIFFocc = 'whole power spectrum in the max modulated voxels or selection (of Nxovels), (pre-post)/pre difference';

readme.Palpha_max       = 'P in the alphamax range in the           max voxel, for averaged alphamax range: Ranksum p for the alphamax 45-90 Hz range';
readme.Palpha_occmax    = 'P in the alphamax range in the occipital max voxel, for averaged alphamax range: Ranksum p for the alphamax 45-90 Hz range';
readme.Palpharange_sel  = 'probability of increase  in the averaged alpha range at max-occ selection of 25 voxels: Ranksum p for the alphamax 45-90 Hz range';
readme.Pow   = 'max power for the max voxel/selection of Nvoxels, calculated from normalized power of type (POSToccmax{con}-PREoccmax{con})./PREoccmax{con}'  
readme.F     = 'freq of the max for the max voxel/selection of Nvoxels' 
readme.WPow  = 'weighted power for the max voxel/selection of Nvoxels, calculated from normalized power of type (POSToccmax{con}-PREoccmax{con})./PREoccmax{con}'
readme.WF    = 'weighted freq for the max voxel/selection of Nvoxels'

readme.selWMAXF = 'the same as selWF, but the selection of N voxels is based on MAX condition'
readme.selWMAXPow = 'the same as selWPow, but the selection of N voxels is based on MAX condition'
readme.Palpharange_sel_MAXselection = 'probability of increase  in the averaged alpha range at max-occ selection of 25 voxelsvased on MAX: Ranksum p for the alphamax 45-90 Hz range';

readme.MAX_ROI = 'region of the maximally modulated voxel: in the whole brain or in occipital selection - occ';
readme.MAX_coord = 'coordinates of N maximal voxels in whole brain/occipital'

readme.Spectrum_pre_and_post = 'full brain spectra, pre and post-stimulus';
readme.DIFF = '(pre-post)/pre at each frequency in selection of Nvoxels max occ voxels'
readme.binsfreqssel = 'bins used wor selWF or maxWF'

readme.NVOXELS_ID = 'id nimbers of 25 voxels around the max voxel';

readme.MAX_VOXEL_OCC_NUMBER ='number of maximally modulated voxel in occipital selection';

readme.PREvisual  =  'mean normalized power in visual areas in prestimulus';
readme.POSTvisual =  'mean normalized power in visual areas in poststimulus';
readme.PREvisual  =  'mean normalized power in calcarines in prestimulus';
readme.POSTvisual =  'mean normalized power in calcarines in poststimulus';

Freq = freq1.freq; % alpha range frequencies

% FDR05occsel, FDR0001occsel,    PREoccsel, POSToccsel, DIFFoccsel,     Palpha_sel,     selPow, selF, selWF, selWPow
% FDR05occmax, FDR0001occmax,     PREoccmax, POSToccmax, DIFFoccmax,     Palpha_max,      maxPow, maxF, maxWF, maxWPow

save (filename, 'readme', 'PREoccsel', 'POSToccsel', 'DIFFoccsel', 'PREoccmax', 'POSToccmax', 'DIFFoccmax', 'Nvoxels', 'FDR05occsel', 'FDR0001occsel', 'FDR05occmax', 'FDR0001occmax',...
      'Freq', 'Palpha_max', 'Palpha_occmax', 'Palpharange_sel',...
      'selPow', 'selF', 'selWF', 'selWPow', 'maxPow', 'maxF', 'maxWF', 'maxWPow', 'MAX_ROI', 'MAX_coord_occ', 'MAX_coord','alphamax', 'SORTED', 'SORTEDocc',...
      'binsfreqssel', 'binsfreqsmax', 'Spectrum_pre', 'Spectrum_post', 'MAX_VOXEL_OCC_NUMBER', 'MAX_VOXEL_NUMBER', 'NVOXELS_ID', ...
      'selWMAXF', 'selWMAXPow', 'Palpharange_sel_MAXselection');
%%
exportToPPTX('saveandclose', PPTXname);



