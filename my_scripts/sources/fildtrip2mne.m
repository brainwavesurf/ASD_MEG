%Save fieldtrip data fo mne python
clear all
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
savepath = '/net/server/data/Archive/aut_gamma/orekhova/KI/Scripts_bkp/Shishkina/KI/Results_Alpha_and_Gamma/';
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
    subj = '0106'; %'SUBJ (s,:); 
    
    savemegto = strcat(savepath, subj);
    epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
    

    % read preprocessed data, 10-17 Hz bandpass
    bp = load(strcat(epofolder, subj, '_preproc_alpha_bp_epochs.mat'));

    % select mag epochs for slow and fast conditions in interstimuli [-0.8 0] and stim [0.4 1.2] period
    cfg = [];
    cfg.channel = 'megmag';
    data_slow = ft_selectdata(cfg, bp.slow_alpha_bp);
    data_fast = ft_selectdata(cfg, bp.fast_alpha_bp); 

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

    epochs_fast1 = data_fast; epochs_fast1.trial = fast{1};
    epochs_fast2 = data_fast; epochs_fast2.trial = fast{2};
    epochs_fast3 = data_fast; epochs_fast3.trial = fast{3};
    epochs_fast4 = data_fast; epochs_fast4.trial = fast{4};
    epochs_fast5 = data_fast; epochs_fast5.trial = fast{5};
    epochs_fast6 = data_fast; epochs_fast6.trial = fast{6};   
    
    epochs_slow1 = data_slow; epochs_slow1.trial = slow{1};
    epochs_slow2 = data_slow; epochs_slow2.trial = slow{2};
    epochs_slow3 = data_slow; epochs_slow3.trial = slow{3};
    epochs_slow4 = data_slow; epochs_slow4.trial = slow{4};    
    epochs_slow5 = data_slow; epochs_slow5.trial = slow{5};
    epochs_slow6 = data_slow; epochs_slow6.trial = slow{6};


    %save freq analysis results
    filename = strcat(savepath, subj, '/', subj, '_fieldtrip_csp.mat');
    
    save(filename, 'epochs_fast1', 'epochs_fast2','epochs_fast3','epochs_fast4','epochs_fast5','epochs_fast6',...
    'epochs_slow1','epochs_slow2', 'epochs_slow3','epochs_slow4','epochs_slow5','epochs_slow6')

end