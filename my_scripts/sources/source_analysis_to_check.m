% take 1 components for fast cond and 6 components for slow cond 
subj = '0106'; 

% load headmodel
headmodel = load('0106_individ_hdm_vol.mat'));
% load leadfield 
leadfield = load('0106_grid_MNI_lf.mat'));
% load mri
mri = load('0106_mri_orig_realigned.mat'));
% load band_passed data
bp = load('0106_preproc_alpha_bp_epochs.mat'));

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
load('0106_csp_analysis.mat'));

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


% for averaging of pow
% pow     = 0;
% ntrial  = 0;

% components for fast cond
for i = 1:6 
    close all
    
    data_comp{i} = data_fast;
    data_comp{i}.trial = fast{i}; % replace MEG data trials with csp (chan x time)
    
    % do timelock and %keep trials 
    cfg                   = [];
    cfg.preproc.demean    = 'yes';    % enable demean to remove mean value from each single trial
    cfg.channel           = 'MEGMAG';
    %cfg.keeptrials        = 'yes';
    MEG_trials_fast{i}         = ft_timelockanalysis(cfg, data_comp{i});
    MEG_trials_fast{i}.cov     = MEG_cov.cov;
    
    
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
    
%     % divide data to single trials 
%     for t = 1:size(MEG_trials{i}.trial,1)
%         
%         MEG_single_trial{i}{t} = MEG_trials{i};
%         
%         % average over trials to make chan_time dimord for every trial
%         cfg = [];
%         cfg.avgoverrpt = 'yes';
%         MEG_single_trial{i}{t} = ft_selectdata(cfg, MEG_single_trial{i}{t});
%         
%         MEG_single_trial{i}{t}.trial = squeeze(MEG_trials{i}.trial(t,:,:));
%         MEG_single_trial{i}{t}.avg = squeeze(MEG_trials{i}.trial(t,:,:));
%         MEG_single_trial{i}{t}.cov = MEG_cov.cov; %add noise cov matrix for each trial
        
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
        cfg.method                  = 'eloreta';                    %specify method 
        cfg.sourcemodel             = leadfield.grid_MNI_lf;        %the precomputed leadfield
        cfg.headmodel               = headmodel.individ_hdm_vol;    %the head model
        %cfg.eloreta.prewhiten       = 'yes';                        %prewhiten data
        %cfg.eloreta.scalesourcecov  = 'yes';                        %scaling the source covariance matrix
        cfg.eloreta.lambda          = 0.05;                         %0.05regularisation parameter - try different values (3)
        cfg.channel                 = 'MEGMAG';
        %source_trials{t}            = ft_sourceanalysis(cfg, MEG_single_trial{i}{t});
        source_trials_fast{i}       = ft_sourceanalysis(cfg, MEG_trials_fast{i});
        %source_changed_fast{i}       = ft_sourceanalysis(cfg, random_MEG_trials_fast{i});   
        
%results for source for each trials       
%       time: [1×401 double]
%        dim: [29 35 30]
%     inside: [30450×1 logical]
%        pos: [30450×3 double]
%     method: 'average'
%        avg: [1×1 struct]
%        cfg: [1×1 struct]

%         %average parameter 'pow' for all trials 
%         source_avg{i}   = source_trials{1}; 
%         pow         = (pow + source_trials{t}.avg.pow)/(ntrial + t);
    end
    
%     % replace pow with average for all trials
%     source_avg{i}.avg.pow = pow;
%     
%     % interpolate data
%     cfg            = [];
%     cfg.parameter  = 'pow';
%     interpolate{i}    = ft_sourceinterpolate(cfg, source_trials{t} , mri.mri_orig_realigned);
% 
%    
%     % plot ortho
%     cfg = [];
%     cfg.method        = 'ortho';
%     cfg.funparameter  = 'pow';
%     ft_sourceplot(cfg,interpolate{i});
%     

%     
%     % save figires
%     saveas(figure(1), [savepath, subj, '/', subj, '_csp_ortho_comp_', num2str(i),'.jpeg']);


% for averaging parameter 'pow' for all trials
% pow     = 0;
% ntrial  = 0;

% components for slow cond
for i = 1:6
    close all
    
    data_comp{i} = data_slow;
    data_comp{i}.trial = slow{i};
    
    % do timelock %for single trials
    cfg                   = [];
    cfg.preproc.demean    = 'yes';    % enable demean to remove mean value from each single trial
    cfg.channel           = 'MEGMAG';
    %cfg.keeptrials        = 'yes';
    MEG_trials_slow{i}         = ft_timelockanalysis(cfg, data_comp{i});
    MEG_trials_slow{i}.cov     = MEG_cov.cov;
    
    
    %check shuffled noise cov
    %MEG_trials_slow{i}.cov = random_MEG_cov;
    
    %shuffle MEG_data
%     random_MEG_trials_slow{i} =  MEG_trials_slow{i};
%     random_MEG_trials_slow{i}.avg = MEG_trials_slow{i}.avg(randperm(size(MEG_trials_slow{i}.avg, 1)), :);

% %%take identical matrix    
% MEG_trials_slow{i}.cov = eye(102); 
  
%     % divide data to single trials 
%     for t = 1:size(MEG_trials{i}.trial,1)
%         
%         MEG_single_trial{i}{t} = MEG_trials{i};
%         MEG_single_trial{i}{t}.avg = MEG_trials{i}.trial(t,:,:);
%         
%         % average over trials to make chan_time dimord for every trial
%         cfg = [];
%         cfg.avgoverrpt = 'yes';
%         MEG_single_trial{i}{t} = ft_selectdata(cfg, MEG_single_trial{i}{t});
%         MEG_single_trial{i}{t}.cov = MEG_cov.cov;  %add noise cov matrix for each trial 
        
    
        % do eloreta for each trial   
        cfg                         = [];
        cfg.method                  = 'eloreta';                    %specify method 
        cfg.sourcemodel             = leadfield.grid_MNI_lf;        %the precomputed leadfield
        cfg.headmodel               = headmodel.individ_hdm_vol;    %the head model
        cfg.eloreta.prewhiten       = 'yes';                        %prewhiten data
        cfg.eloreta.scalesourcecov  = 'yes';                        %scaling the source covariance matrix
        cfg.eloreta.lambda          = 0.05;                         %regularisation parameter - try different values (3)
        cfg.channel                 = 'MEGMAG';
        %source_trials{t}            = ft_sourceanalysis(cfg, MEG_single_trial{i}{t});
        source_trials_slow{i}       = ft_sourceanalysis(cfg, MEG_trials_slow{i});
        %source_changed_slow{i}       = ft_sourceanalysis(cfg, random_MEG_trials_slow{i});   
        

%         source_avg{i} = source_trials{1};
%         pow       = (pow + source_trials{t}.avg.pow)/(ntrial + t);
%     
end

%localisation of difference
for c = 1:6
% replace pow with difference between power in two conditions
    source_diff{c} = source_trials_fast{c};
    %source_diff{c}.avg.pow = source_changed_fast{1}.avg.pow - source_changed_slow{1}.avg.pow;
    source_diff{c}.avg.pow = source_trials_fast{c}.avg.pow - source_trials_slow{c}.avg.pow;

    % interpolate data
    cfg            = [];
    cfg.parameter  = 'pow';
    interpolate{c}    = ft_sourceinterpolate(cfg, source_diff{c}, mri.mri_orig_realigned);
    %interpolate{i}    = ft_sourceinterpolate(cfg, source_diff{c}, mri.mri_orig_realigned);

    % plot ortho
    cfg = [];
    cfg.method        = 'ortho';
    cfg.funparameter  = 'pow';
    ft_sourceplot(cfg,interpolate{c});
    
    %saveas(figure(2), [savepath, subj, '/', subj, '_csp_ortho_scp1fast_csp6slow.jpeg']);
end
    
