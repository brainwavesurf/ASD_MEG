%Source localization of MEG data after CSP by 'eloreta' for one subject
clear all
close all
% path info
fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(fieldtripfolder, path);
ft_defaults;

path('/home/a_shishkina/externals/files', path);
path('/home/a_shishkina/data/KI/SUBJECTS/0106/ICA_nonotch_crop', path);
path('/home/a_shishkina/fieldtrip/template/anatomy', path)

realdatapath = '/home/a_shishkina/data/KI/SUBJECTS/';
data_path = '/home/a_shishkina/data/KI/FT_beamformer/';
mrifolder = strcat( '/mri_linwarp_', num2str(6), 'mm_', 'brthr',  '0.9/');
savepath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';
%%

% load subj info
     
SUBJ_NT = [ '0101'; '0102'; '0103'; '0104'; '0105'; '0135'; '0136';...  
            '0137'; '0138'; '0140'; '0158'; '0162'; '0163'; '0178';...
            '0179'; '0255'; '0257'; '0348'; '0378'; '0379'; '0384']; 
        
SUBJ_ASD = ['0106'; '0107'; '0139'; '0141'; '0159'; '0160'; '0161';...  
            '0164'; '0253'; '0254'; '0256'; '0273'; '0274'; '0275';...
            '0276'; '0346'; '0347'; '0349'; '0351'; '0358';...
            '0380'; '0381'; '0382'; '0383'];  
%without '0357';
SUBJ = [SUBJ_ASD; SUBJ_NT];

%% loop for all subjects
for s=1: size (SUBJ,1)
    
    close all
    subj = '0106';%SUBJ (s,:); 
    
    savemegto = strcat(savepath, subj);
    epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');

    % load headmodel
    headmodel = load(strcat(data_path, subj, mrifolder, subj, '_individ_hdm_vol.mat'));
    % load leadfield 
    leadfield = load(strcat(data_path, subj, mrifolder, subj, '_grid_MNI_lf.mat'));
    % load mri
    mri = load(strcat(data_path, subj, mrifolder, subj, '_mri_orig_realigned.mat')); 
    % 
    % read preprocessed data, 10-14 Hz bandpass
    epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
    bp = load(strcat(epofolder, subj, '_preproc_alpha_bp_epochs.mat'));
    
    % select mag epochs for slow and fast conditions in interstimuli [-0.8 0] and stim [0.4 1.2] period
    cfg = [];
    cfg.channel = 'MEGMAG';
    data_slow = ft_selectdata(cfg, bp.slow_alpha_bp);
    data_fast = ft_selectdata(cfg, bp.fast_alpha_bp);

    data_slow_post = ft_selectdata(cfg, bp.slow_alpha_post);
    data_fast_post = ft_selectdata(cfg, bp.fast_alpha_post);

    % use concatenated data for noise covariance matrix calculation, as alpha is suppressed here
    data_for_cov = ft_appenddata([],data_slow_post, data_fast_post);

    % load csp data for two conditions
    load(strcat(savepath, subj, '/', subj, '_csp_analysis.mat'));

    % converted the Xcsp_fast and Xcsp_slow to MEG time series
    for s = 1:6
        A_mat = A1;
        for i = 1:6
            if i~=s
                A_mat(i,:) = 0;
            end
            A{s} = A_mat; % leave only one component, the rest ones are 0s 
        end
        %for fast
        trial_fast = zeros(size(Xcsp_fast,1), size(Xcsp_fast,2), 102);
        for j = 1:size(Xcsp_fast,1)
            trial_fast(j,:,:) = squeeze(Xcsp_fast(j,:,:))*A{s};
            csp_data_fast{j} = squeeze(trial_fast(j,:,:))';
        end
        fast{s} = csp_data_fast; % first three components for fast cond
        %for slow
        trial_slow = zeros(size(Xcsp_slow,1), size(Xcsp_slow,2), 102);
        for j = 1:size(Xcsp_slow,1)
            trial_slow(j,:,:) = squeeze(Xcsp_slow(j,:,:))*A{s};
            csp_data_slow{j} = squeeze(trial_slow(j,:,:))';
        end
        slow{s} = csp_data_slow; % last three components for slow cond
    end
end
% timelock for noise covariance matrix calculation on MEG data for two condition in stim interval 
cfg                   = [];
cfg.preproc.demean    = 'yes';      % enable demean to remove mean value from each single trial
cfg.channel           = 'MEGMAG';
cfg.covariance        = 'yes';      % calculate covariance matrix of the data
cfg.covariancewindow  = [0.4 1.2];  % calculate the covariance matrix for a specific time window
MEG_cov               = ft_timelockanalysis(cfg, data_for_cov);


%for averaging of pow
pow     = 0;
ntrial  = 0;

% components for fast cond
for i = 1:2
    close all
    
    data_comp{i} = data_fast;
    data_comp{i}.trial = fast{i}; % replace MEG data trials with csp (chan x time)
    
    % do timelock and %keep trials 
    cfg                   = [];
    cfg.preproc.demean    = 'yes';    % enable demean to remove mean value from each single trial
    cfg.channel           = 'MEGMAG';
    cfg.keeptrials        = 'yes';
    MEG_trials_fast{i}         = ft_timelockanalysis(cfg, data_comp{i});
    %MEG_trials_fast{i}.cov     = MEG_cov.cov;
    
    
%     %shuffle cov
%     D = diag(MEG_cov.cov);	
%     random_MEG_cov = MEG_cov.cov(randperm(size(MEG_cov.cov, 1)), :);    
%     random_MEG_cov  = random_MEG_cov  - diag(diag(random_MEG_cov )) + diag(D);

    %%take identical matrix
    %MEG_trials_fast{i}.cov = eye(102); 
    %check shuffled noise cov
    %MEG_trials_fast{i}.cov = random_MEG_cov;
    
    %shuffle MEG_data
    %random_MEG_trials_fast{i} =  MEG_trials_fast{i};
    %random_MEG_trials_fast{i}.avg = MEG_trials_fast{i}.avg(randperm(size(MEG_trials_fast{i}.avg, 1)), :);
    
    % divide data to single trials 
    for t = 1:2%size(MEG_trials{i}.trial,1)
        
        MEG_single_trial{i}{t} = MEG_trials_fast{i};
        
        % average over trials to make chan_time dimord for every trial
        cfg = [];
        cfg.avgoverrpt = 'yes';
        MEG_single_trial{i}{t} = ft_selectdata(cfg, MEG_single_trial{i}{t});
        
        MEG_single_trial{i}{t}.trial = squeeze(MEG_trials_fast{i}.trial(t,:,:));
        MEG_single_trial{i}{t}.avg = squeeze(MEG_trials_fast{i}.trial(t,:,:));
        MEG_single_trial{i}{t}.cov = MEG_cov.cov; %add noise cov matrix for each trial
        
%         plot(MEG_trials{i}.avg(1,:), '-r')
%         hold on
%         plot(MEG_single_trial{i}{t}.avg(1,:))

%results after timelock for each trials and avg after
%       time: [1×401 double]
%      label: {102×1 cell}
%       grad: [1×1 struct]
%      trial: [102×401 double]
%        cfg: [1×1 struct]
%     dimord: 'chan_time'
%        avg: [102×401 double]
%        cov: [102×102 double]
        
        % do eloreta for each trial
        cfg                         = [];
        cfg.method                  = 'mne';                    %specify method 
        cfg.grad = MEG_single_trial{i}{t}.grad;   % my addition
        cfg.channel = {'MEGMAG'}; % my addition
        cfg.sourcemodel             = leadfield.grid_MNI_lf;        %the precomputed leadfield
        cfg.headmodel               = headmodel.individ_hdm_vol;    %the head model
        cfg.mne.prewhiten       = 'yes';                        %prewhiten data
        cfg.mne.scalesourcecov  = 'yes';                        %scaling the source covariance matrix
        cfg.mne.lambda          = 0.05;                         %0.05regularisation parameter - try different values (3)
        cfg.channel                 = 'MEGMAG';
        source_trials{t}            = ft_sourceanalysis(cfg, MEG_single_trial{i}{t});
        %source_trials_fast{i}       = ft_sourceanalysis(cfg, MEG_trials_fast{i});
        %source_changed_fast{i}       = ft_sourceanalysis(cfg, random_MEG_trials_fast{i});   
        
%results for source for each trials       
%       time: [1×401 double]
%        dim: [29 35 30]
%     inside: [30450×1 logical]
%        pos: [30450×3 double]
%     method: 'average'
%        avg: [1×1 struct]
%        cfg: [1×1 struct]

        %average parameter 'pow' for all trials 
        source_avg{i}   = source_trials{1}; 
        pow         = (pow + source_trials{t}.avg.pow)/(ntrial + t);
    end
    
    % replace pow with average for all trials
    source_avg{i}.avg.pow = pow;
%     
    % interpolate data
    cfg            = [];
    cfg.parameter  = 'pow';
    interpolate{i}    = ft_sourceinterpolate(cfg,  source_avg{i} , mri.mri_orig_realigned);

%    
    % plot ortho
    cfg = [];
    cfg.method        = 'ortho';
    cfg.funparameter  = 'pow';
    ft_sourceplot(cfg,interpolate{i});
%     
%     % spatially normalize the anatomy and functional data to MNI coordinates
%     cfg = [];
%     cfg.nonlinear = 'no';
%     normalize{i} = ft_volumenormalise(cfg, interpolate{i});
%     
%     % surface plot
%     cfg = [];
%     cfg.method         = 'surface';
%     cfg.funparameter   = 'pow';
%     cfg.maskparameter  = cfg.funparameter;
%     %cfg.funcolorlim    = [0.0 2.4e-16];
%     cfg.funcolormap    = 'jet';
%     %cfg.opacitylim     = [0.0 2.4e-16];
%     cfg.opacitymap     = 'rampup';
%     cfg.projmethod     = 'nearest';
%     cfg.surffile       = 'surface_white_both.mat'; % Cortical sheet from canonical MNI brain
%     ft_sourceplot(cfg, normalize{i});
%     view ([45 30])
%     
%     % save figires
%     saveas(figure(1), [savepath, subj, '/', subj, '_csp_ortho_comp_', num2str(i),'.jpeg']);
%     saveas(figure(2), [savepath, subj, '/', subj, '_csp_surface_comp_', num2str(i),'.jpeg']);


% for averaging parameter 'pow' for all trials
pow     = 0;
ntrial  = 0;

% components for slow cond
for i = 1:6
    close all
    
    data_comp{i} = data_slow;
    data_comp{i}.trial = slow{i};
    
    % do timelock %for single trials
    cfg                   = [];
    cfg.preproc.demean    = 'yes';    % enable demean to remove mean value from each single trial
    cfg.channel           = 'MEGMAG';
    cfg.keeptrials        = 'yes';
    MEG_trials_slow{i}         = ft_timelockanalysis(cfg, data_comp{i});
    %MEG_trials_slow{i}.cov     = MEG_cov.cov;
    
    
    %check shuffled noise cov
    %MEG_trials_slow{i}.cov = random_MEG_cov;
    
    %shuffle MEG_data
%     random_MEG_trials_slow{i} =  MEG_trials_slow{i};
%     random_MEG_trials_slow{i}.avg = MEG_trials_slow{i}.avg(randperm(size(MEG_trials_slow{i}.avg, 1)), :);
% %%take identical matrix    
% MEG_trials_slow{i}.cov = eye(102); 
  
    % divide data to single trials 
    for t = 1:2 %size(MEG_trials{i}.trial,1)
        
        MEG_single_trial{i}{t} = MEG_trials_slow{i};
        MEG_single_trial{i}{t}.avg = MEG_trials_slow{i}.trial(t,:,:);
        
        % average over trials to make chan_time dimord for every trial
        cfg = [];
        cfg.avgoverrpt = 'yes';
        MEG_single_trial{i}{t} = ft_selectdata(cfg, MEG_single_trial{i}{t});
        MEG_single_trial{i}{t}.cov = MEG_cov.cov;  %add noise cov matrix for each trial 
        
        %  % do eloreta for each trial
        cfg                         = [];
        cfg.method                  = 'mne';                    %specify method 
        cfg.grad = MEG_single_trial{i}{t}.grad;   % my addition
        cfg.channel = {'MEGMAG'}; % my addition
        cfg.sourcemodel             = leadfield.grid_MNI_lf;        %the precomputed leadfield
        cfg.headmodel               = headmodel.individ_hdm_vol;    %the head model
        cfg.mne.prewhiten           = 'yes';                        %prewhiten data
        cfg.mne.scalesourcecov      = 'yes';                        %scaling the source covariance matrix
        cfg.mne.lambda              = 0.05;                         %0.05regularisation parameter - try different values (3)
        cfg.channel                 = 'MEGMAG';
        source_trials_slow{t}       = ft_sourceanalysis(cfg, MEG_single_trial{i}{t});
        
        % do eloreta for each trial   
        cfg                         = [];
        cfg.method                  = 'eloreta';                    %specify method 
        cfg.sourcemodel             = leadfield.grid_MNI_lf;        %the precomputed leadfield
        cfg.headmodel               = headmodel.individ_hdm_vol;    %the head model
        %cfg.eloreta.prewhiten       = 'yes';                        %prewhiten data
        %cfg.eloreta.scalesourcecov  = 'yes';                        %scaling the source covariance matrix
        cfg.eloreta.lambda          = 0.05;                         %regularisation parameter - try different values (3)
        cfg.channel                 = 'MEGMAG';
        source_trials_slow{t}            = ft_sourceanalysis(cfg, MEG_single_trial{i}{t});
        %source_trials_slow{i}       = ft_sourceanalysis(cfg, MEG_trials_slow{i});
        %source_changed_slow{i}       = ft_sourceanalysis(cfg, random_MEG_trials_slow{i});   
        

        source_avg_slow{i} = source_trials_slow{1};
        pow       = (pow + source_trials_slow{t}.avg.pow)/(ntrial + t);
    
    end

% replace pow with average for all trials
    source_avg_slow{i}.avg.pow = pow;
%     
    % interpolate data
    cfg            = [];
    cfg.parameter  = 'pow';
    interpolate_slow{i}    = ft_sourceinterpolate(cfg,  source_avg_slow{i} , mri.mri_orig_realigned);

%    
    % plot ortho
    cfg = [];
    cfg.method        = 'ortho';
    cfg.funparameter  = 'pow';
    ft_sourceplot(cfg,interpolate_slow{i});
    
    
for c = 1:6
% replace pow with difference between power in two conditions
    source_diff{c} = source_trials_fast{1};
    %source_diff{c}.avg.pow = source_changed_fast{1}.avg.pow - source_changed_slow{1}.avg.pow;
    source_diff{c}.avg.pow = source_trials{t}.avg.pow - source_trials_slow{t}.avg.pow;
    
    source_trials{t}.avg.pow = nanmean(source_trials{t}.avg.pow,2);
    % interpolate data
    cfg            = [];
    cfg.parameter  = 'pow';
    interpolate{i}    = ft_sourceinterpolate(cfg, source_trials{t}, mri.mri_orig_realigned);
    %interpolate{i}    = ft_sourceinterpolate(cfg, source_diff{c}, mri.mri_orig_realigned);

    % plot ortho
    cfg = [];
    cfg.method        = 'ortho';
    cfg.funparameter  = 'pow';
    ft_sourceplot(cfg,interpolate{i});
    
    
    %saveas(figure(2), [savepath, subj, '/', subj, '_csp_ortho_scp1fast_csp6slow.jpeg']);
end
    
%     % spatially normalize the anatomy and functional data to MNI coordinates
%     cfg = [];
%     cfg.nonlinear = 'no';
%     normalize{i} = ft_volumenormalise(cfg, interpolate{i});
    
%     % surface plot
%     cfg = [];
%     cfg.method         = 'surface';
%     cfg.funparameter   = 'pow';
%     cfg.maskparameter  = cfg.funparameter;
%     %cfg.funcolorlim    = [0.0 2.4e-16];
%     cfg.funcolormap    = 'jet';
%     %cfg.opacitylim     = [0.0 2.4e-16];
%     cfg.opacitymap     = 'rampup';
%     cfg.projmethod     = 'nearest';
%     cfg.surffile       = 'surface_white_both.mat'; % Cortical sheet from canonical MNI brain
%     ft_sourceplot(cfg, normalize{i});
%     view ([45 30])
%     
%     %save figires
     %saveas(figure(2), [savepath, subj, '/', subj, '_csp_ortho_scp1fast_csp6slow', num2str(i),'.jpeg']);
%     saveas(figure(2), [savepath, subj, '/', subj, '_csp_surface_comp_', num2str(i),'.jpeg']);
% end
% 
% %compare scale
% figure(2)
% subplot_tight(1,2,1,[0.08 0.08])
% h1 = heatmap(MEG_single_trial{i}{t}.cov)
% h1.Title = 'Noise cov';
% h1.XLabel = 'channel';
% h1.YLabel = 'channel';
% subplot_tight(1,2,2,[0.08 0.08])
% h2 = heatmap(MEG_single_trial{i}{t}.avg(:,200:401))
% h2.Title = 'MEG data';
% h2.YLabel = 'channel';
% h2.XLabel = 'time points';
% saveas(figure(1), [savepath, subj, '/', subj, '_noisecov_vs_megdata.jpeg']);