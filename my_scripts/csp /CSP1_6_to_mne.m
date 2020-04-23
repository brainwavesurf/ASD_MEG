%Save fieldtrip data fo mne python
close all;
% path info
fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(fieldtripfolder, path);
ft_defaults;

path('/home/a_shishkina/externals/files', path);
path('/home/a_shishkina/fieldtrip/template/anatomy', path)

realdatapath = '/net/server/data/Archive/aut_gamma/orekhova/KI/SUBJECTS/';
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

SUBJ = [SUBJ_ASD; SUBJ_NT];

%% loop for all subjects
for s=1: size(SUBJ,1)
    
    close all
    subj = SUBJ(s,:); 

    % read preprocessed data, 10-17 Hz bandpass
    load(strcat(savepath, subj, '/', subj, '_preproc_alpha_10_17_epochs.mat'));

    % select mag epochs for slow and fast conditions in interstimuli [-0.8 0] and stim [0.4 1.2] period
    cfg = [];
    cfg.channel = 'megmag';
    data_slow = ft_selectdata(cfg, slow_alpha_isi);
    data_fast = ft_selectdata(cfg, fast_alpha_isi); 
    % load csp data for two conditions
    load(strcat(savepath, subj, '/', subj, '_csp_analysis_1_6.mat'));

    % converted the Xcsp_fast and Xcsp_slow to MEG time series
    A = cell(1,12);
    csp_data_fast = cell(1,size(Xcsp_fast,1));
    fast = cell(1,12);
    csp_data_slow = cell(1,size(Xcsp_slow,1));
    slow = cell(1,12);
    epochs_fast = cell(1,12);
    epochs_slow = cell(1,12);
    for n = 1:12
        A_mat = A1;
        for i = 1:12
            if i~=n
                A_mat(i,:) = 0;
            end
            A{n} = A_mat; % leave only one component, the rest ones are 0s 
        end
        %for fast
        trial_fast = zeros(size(Xcsp_fast,1), size(Xcsp_fast,2), 102);
        for j = 1:size(Xcsp_fast,1)
            trial_fast(j,:,:) = squeeze(Xcsp_fast(j,:,:))*A{n};
            csp_data_fast{j} = squeeze(trial_fast(j,:,:))';
        end
        fast{n} = csp_data_fast; % first three components for fast cond
        %for slow
        trial_slow = zeros(size(Xcsp_slow,1), size(Xcsp_slow,2), 102);
        for j = 1:size(Xcsp_slow,1)
            trial_slow(j,:,:) = squeeze(Xcsp_slow(j,:,:))*A{n};
            csp_data_slow{j} = squeeze(trial_slow(j,:,:))';
        end
        slow{n} = csp_data_slow; % last three components for slow cond
        
        epochs_fast{n} = data_fast; epochs_fast{n}.trial = fast{n};
        epochs_slow{n} = data_slow; epochs_slow{n}.trial = slow{n};
    end

    %save freq analysis results
    filename = strcat(savepath, subj, '/', subj, '_fieldtrip_csp_1_6.mat');
    
    save(filename, 'epochs_fast', 'epochs_slow')
end

%plot time series
for i=7:12
    figure(1)
    subplot(2,3,i-6)
    plot(fast{i}{1}(1,1:300), '-r')
    hold on
    plot(slow{i}{1}(1,1:300), '-b')
    legend('fast', 'slow');
    title(strcat('csp', num2str(90+i)))
end
saveas(figure(1), [savepath, '1_results/CSP_matlab_plots/', subj, '_csp_97_102_timeseries.jpeg']);  

%double check to estimate the spectral power for all eigenvalues
file = [epochs_fast, epochs_slow];
  
for num = 7:12
    cfg = [];
    cfg.method       = 'mtmfft';
    cfg.output       = 'pow'; 
    cfg.taper        = 'hanning'; 
    cfg.pad          = 10; 
    cfg.foilim       = [6 20];
    cfg.tapsmofrq    = 3; 
    fft_fast  = ft_freqanalysis(cfg, epochs_fast{num}); 
    fft_slow   = ft_freqanalysis(cfg, epochs_slow{num});

    cfg = [];
    cfg.avgoverchan = 'yes';
    fft_fast = ft_selectdata(cfg,fft_fast);
    fft_slow = ft_selectdata(cfg,fft_slow);

    figure(1)
    subplot(2,3,num-6)
    plot(fft_fast.freq, fft_fast.powspctrm, '-r')
    hold on;
    plot(fft_slow.freq, fft_slow.powspctrm, '-b')
    legend('fast', 'slow'); 
    title(strcat('csp', num2str(90+num)))
end
saveas(figure(1), [savepath, '1_results/CSP_matlab_plots/', subj, '_spectral_power_csp_97_102_timeseries.jpeg']);