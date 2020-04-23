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
    load(strcat(savepath, subj, '/', subj, '_csp_analysis_4_6.mat'));

    % converted the Xcsp_fast and Xcsp_slow to MEG time series
    A = cell(1,3);
    csp_data_fast = cell(1,size(Xcsp_fast,1));
    fast = cell(1,3);
    csp_data_slow = cell(1,size(Xcsp_slow,1));
    slow = cell(1,3);
    for n = 1:3
        A_mat = A1;
        for i = 1:3
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
    end

    epochs_fast1 = data_fast; epochs_fast1.trial = fast{1};
    epochs_fast2 = data_fast; epochs_fast2.trial = fast{2};
    epochs_fast3 = data_fast; epochs_fast3.trial = fast{3}; 
    
    epochs_slow1 = data_slow; epochs_slow1.trial = slow{1};
    epochs_slow2 = data_slow; epochs_slow2.trial = slow{2};
    epochs_slow3 = data_slow; epochs_slow3.trial = slow{3};

    %save freq analysis results
    filename = strcat(savepath, subj, '/', subj, '_fieldtrip_csp_4_6.mat');
    
    save(filename, 'epochs_fast1', 'epochs_fast2','epochs_fast3',...
                   'epochs_slow1','epochs_slow2', 'epochs_slow3')

end

for i=[1,2,3]
    figure(1)
    subplot(1,3,i)
    plot(fast{i}{1}(1,1:300), '-r')
    hold on
    plot(slow{i}{1}(1,1:300), '-b')
    legend('fast', 'slow');
    title(strcat('csp', num2str(i+3)))
end
saveas(figure(1), [savepath, '1_results/CSP_matlab_plots/', subj, '_csp_4_6_timeseries.jpeg']);  

%double check to estimate the spectral power for all eigenvalues
file = [epochs_fast1, epochs_fast2, epochs_fast3,...
        epochs_slow1, epochs_slow2, epochs_slow3];
    
for num = 1:6
    cfg = [];
    cfg.method       = 'mtmfft';
    cfg.output       = 'pow'; 
    cfg.taper        = 'hanning'; 
    cfg.pad          = 10; 
    cfg.foilim       = [10 17];
    cfg.tapsmofrq    = 3; 
    fft_fast  = ft_freqanalysis(cfg, file(num)); 
    fft_slow   = ft_freqanalysis(cfg, file(num));

    cfg = [];
    cfg.avgoverchan = 'yes';
    fft_fast = ft_selectdata(cfg,fft_fast);
    fft_slow = ft_selectdata(cfg,fft_slow);


    figure(1)
    subplot(2,3,num)
    plot(fft_fast.freq, fft_fast.powspctrm, '-r')
    hold on;
    plot(fft_slow.freq, fft_slow.powspctrm, '-b')
    legend('fast', 'slow'); 
    title(strcat('csp', num2str(num)))
end
saveas(figure(1), [savepath, subj, '/', subj, '_spectral_power_csp_timeseries.jpeg']);