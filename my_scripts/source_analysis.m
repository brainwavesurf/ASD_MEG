%Source localization of MEG data after CSP by 'eloreta' for one subject

% path info
fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(fieldtripfolder, path);
ft_defaults;

path('/home/a_shishkina/externals/files', path);
path('/home/a_shishkina/data/KI/SUBJECTS/0106/ICA_nonotch_crop', path);

realdatapath = '/home/a_shishkina/data/KI/SUBJECTS/';
data_path = '/home/a_shishkina/data/KI/FT_beamformer/';
mrifolder = strcat( '/mri_linwarp_', num2str(6), 'mm_', 'brthr',  '0.9/');
savepath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';

% for one subj
subj = '0106';

% load headmodel
headmodel = load(strcat(data_path, subj, mrifolder, subj, '_individ_hdm_vol.mat'));
% load leadfield 
leadfield = load(strcat(data_path, subj, mrifolder, subj, '_grid_MNI_lf.mat'));
% load mri
mri = load(strcat(data_path, subj, mrifolder, subj, '_mri_orig_realigned.mat')); 

% read preprocessed data
epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
bp = load(strcat(epofolder, subj, '_preproc_alpha_bp_epochs.mat'));

% select mag epochs separately for slow and fast conditions in interstimuli and stim period
cfg = [];
cfg.channel = 'MEGMAG';
data_slow = ft_selectdata(cfg, bp.slow_alpha_bp);
data_fast = ft_selectdata(cfg, bp.fast_alpha_bp);

% use this data for noise covariance matrix calculation
data_slow_post = ft_selectdata(cfg, bp.slow_alpha_post);
data_fast_post = ft_selectdata(cfg, bp.fast_alpha_post);
data_for_cov = ft_appenddata([],data_slow_post, data_fast_post);

% plot one trial
figure;plot(data_segmented_alpha.trial{10}')

% load csp data
load(strcat(savepath, subj, '/', subj, '_csp_analysis.mat'));

% converted the Xcsp_fast and Xcsp_slow to MEG time series
for s = 1:6
    if s<=3
        A_mat = A1;
        for i = 1:6
            if i~=s
                A_mat(i,:) = 0;
            end
            A{s} = A_mat;
        end

        for j = 1:size(Xcsp_fast,1)
            trial_fast(j,:,:) = squeeze(Xcsp_fast(j,:,:))*A{s};
            csp_data_fast{j} = squeeze(trial_fast(j,:,:))';
        end
        fast{s} = csp_data_fast; % first three components for fast cond
    end
    
    else
        A_mat = A1;
        for i = 1:6
            if i~=s
                A_mat(i,:) = 0;
            end
            A{s} = A_mat;
        end

        for j = 1:size(csp.Xcsp_slow,1)
            trial_slow(j,:,:) = squeeze(csp.Xcsp_slow(j,:,:))*A{s};
            csp_data_slow{j} = squeeze(trial_slow(j,:,:))';
        end
        slow{s} = csp_data_slow; % last three components for slow cond
end

% timelock for noise covariance matrix calculation on MEG data for two condition in stim interval 
cfg                   = [];
cfg.preproc.demean    = 'yes';      % enable demean to remove mean value from each single trial
cfg.channel           = 'MEGMAG';
cfg.covariance        = 'yes';      % calculate covariance matrix of the data
cfg.covariancewindow  = [0.4 1.2];  % calculate the covariance matrix for a specific time window
MEG_cov               = ft_timelockanalysis(cfg, data_for_cov);

% replace trials with CSP components, the first three comp for fast cond, the last three for slow cond
for i = 1:3
    data_comp{i} = data_fast;
    data_comp{i}.trial = fast{i};
    
    %do timelock for single trials
    cfg                   = [];
    cfg.preproc.demean    = 'yes';    % enable demean to remove mean value from each single trial
    cfg.channel           = 'MEGMAG';
    cfg.keeptrials        = 'yes';
    MEG_trials{i}         = ft_timelockanalysis(cfg, data_comp{i});
    
    %create noise cov matrix for keeptrials case, same for each trial
    for j = 1:size(MEG_trials{1}.trial,1)
        noise_cov(j,:,:)  = MEG_cov.cov;
    end
    MEG_trials{i}.cov = noise_cov;
    
    %divide data to single trials 
    MEG_single_trial = MEG_trials{i};
    for t = 1:size(MEG_trials{1}.trial,1)
        MEG_single_trial.trial = MEG_trials{i}.trial(j,:,:);
        MEG_single_trl{t} = MEG_single_trial;
    end
    data_single_trial{i} = MEG_single_trl;
end

%do eloreta for each trial 
for t = 1:size(MEG_trials{1}.trial,1)
    cfg                         = [];
    cfg.method                  = 'eloreta';                    %specify method 
    cfg.sourcemodel             = leadfield.grid_MNI_lf;        %the precomputed leadfield
    cfg.headmodel               = headmodel.individ_hdm_vol;    %the head model
    cfg.eloreta.prewhiten       = 'yes';                        %prewhiten data
    cfg.eloreta.scalesourcecov  = 'yes';                        %scaling the source covariance matrix
    cfg.eloreta.lambda          = 0.05;                         %regularisation parameter - try different values (3)
    cfg.channel                 = 'MEGMAG';
    source_trials{t}            = ft_sourceanalysis(cfg, data_single_trial{1}{t});

    %average parameter 'pow' for all trials
    pow     = 0;
    ntrial  = 0;
    source_avg          = source_trials{1}; 
    source_avg.avg.pow  = (pow + source_trials{t}.avg.pow)/(ntrial + t);
end

%interpolate data
cfg            = [];
cfg.parameter  = 'pow';
interpolate    = ft_sourceinterpolate(cfg, source_avg, mri.mri_orig_realigned);

%plot
cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'pow';
ft_sourceplot(cfg,interpolate);

saveas(figure(1), [savepath, subj, '/', subj, '_csp_ortho_comp1.jpeg']);
%%

% for i = 1:2
% trial(i,:,:) = source_trials.trial(i).pow;
% end
% s1 = source_slow;
% s1.trial = trial;
% 
% cfg = [];
% cfg.trials = 1:2;
% cfg.avgoverrpt = 'yes';
% source_slow_avg = ft_selectdata(cfg,s1);
% source_fast_avg = ft_selectdata(cfg,source_fast);
% 
% 
% figure;plot(MEG_fast.time, MEG_fast.avg(1,:))
% hold on
% plot(MEG_fast.time, MEG_slow.avg(1,:))
% 
% figure;plot(source_slow.avg.pow) % with MEG_avg.cov
% hold on
% plot(source_fast.avg.pow) % with MEG_avg_alpha.cov
% 
% cfg            = [];
% cfg.parameter  = 'pow';
% interpolate_slow = ft_sourceinterpolate(cfg, source_slow_avg, mri.mri_orig_realigned);
% interpolate_fast = ft_sourceinterpolate(cfg, source_fast, mri.mri_orig_realigned);
% 
% cfg = [];
% cfg.method        = 'ortho';
% cfg.funparameter  = 'pow';
% ft_sourceplot(cfg,interpolate_slow);
% ft_sourceplot(cfg,interpolate_fast);
% 
% saveas(figure(1), [savepath, subj, '/', subj, '_csp_ortho_slow.jpeg']);
% saveas(figure(2), [savepath, subj, '/', subj, '_csp_ortho_fast.jpeg']);
% 
% % spatially normalize the anatomy and functional data to MNI coordinates
% cfg = [];
% cfg.nonlinear = 'no';
% source_norm_slow = ft_volumenormalise(cfg, interpolate_slow);
% source_norm_fast = ft_volumenormalise(cfg, interpolate_fast);
% 
% % plot multiple 2D axial slices
% cfg = [];
% cfg.method        = 'slice';
% cfg.funparameter  = 'pow';
% %cfg.funcolorlim   = [0.0 1.5e-13];
% ft_sourceplot(cfg, source_norm_slow);
% ft_sourceplot(cfg, source_norm_fast);
% 
% saveas(figure(3), [savepath, subj, '/', subj, 'csp_slice_slow.jpeg']);
% saveas(figure(4), [savepath, subj, '/', subj, 'csp_slice_fast.jpeg']);
% 
% 
% 
% 
