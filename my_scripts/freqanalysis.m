% Load data epochs 
% Load info about the preceding trials

% Calculate time-freq representations by 
% 1) mtmconvol with Hanning taper with varying length of window
% 2) wavelet method with enhancement of the sensitivity to the non-phase locked activity
% 3) mtmfft with multitapers 
% 4) mtmfft with Hanning tapper

% Plot the time-frequency representation (TRF) for grad and mag in the each condition for
% 1) multiple channels (only grad)
% 2) for one channel with max power according to the previous plot (2043)
% 3) time-frequency topoplot for alpha range in interstimulus [-0.8 0] 

% Calculate and plot Fast-Slow power differences for grad and mag

%%
clear all;
close all;

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(fieldtripfolder, path);
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);
realdatapath = '/home/a_shishkina/data/KI/SUBJECTS/';
savepath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';

%%
%add list of subjects:
%SUBJ = ['0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0136'; '0137'; '0138'; '0139'; '0140'; '0141'; '0158'; '0159'; '0160'; '0161'; '0162'; '0163'; '0164'; '0178'; '0179'; '0253'; '0254'; '0255'; '0256'; '0257'; '0259'; '0273'; '0274'; '0275'; '0276'; '0277'; '0346'; '0347'; '0348'; '0350'; '0351'; '0357'; '0358'; '0378'; '0380'; '0381'; '0382'; '0383'; '0384'; '0385']; 
SUBJ = ['0101'];
%%
for s=1: size (SUBJ,1)
    close all
    subj = SUBJ (s,:); 
    savemegto = strcat(savepath, subj);
    epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
    
    %% Load unfiltered epochs. If you need to filter the data (e.g. for LCMV), import raw, not epochs.
    ep_fiff_file = strcat(epofolder, subj, '-noerror-lagcorrected-epo.fif')
    hdr = ft_read_header(ep_fiff_file);
    
    cfg = [];  
    cfg.dataset = ep_fiff_file;
    cfg.channel={'meg'};
    epochs = ft_preprocessing(cfg);
    
    %% Calculating power in sensors, find max weighted freq, freq range and P of alpha enhancement in this freq range.
    
    % posterior sensors, you may use them later for extracting alpha peak power
    %post_sens = {'MEG1932',  'MEG1922', 'MEG2042',  'MEG2032',  'MEG2112', 'MEG2122',  'MEG2342', 'MEG2332',  'MEG1732', 'MEG1942', 'MEG1912', 'MEG2012', 'MEG2022', 'MEG2312', 'MEG2322', 'MEG2512',...
    %     'MEG1933',  'MEG1923', 'MEG2043',  'MEG2033',  'MEG2113', 'MEG2123',  'MEG2343', 'MEG2333',  'MEG1733', 'MEG1943', 'MEG1913', 'MEG2013', 'MEG2023', 'MEG2313', 'MEG2323', 'MEG2513'};
    %% load group preceding events in order to select the epochs later on
    load ([ savemegto, '/', subj, '_info.mat'])
    ev1 = find(allinfo.prev_stim_type==2); %epochs following slow
    ev2 = find(allinfo.prev_stim_type==4); %medium
    ev3 = find(allinfo.prev_stim_type==8); %fast
    prev = {ev1, ev2, ev3};
    
    %%    
    for con=1:3 % for conditions
        %select epochs according to preceding trials
        cfg = [];
        cfg.trials = prev{con}; 
        dataCon = ft_selectdata(cfg, epochs);
        
        %select interstimulus interval
        cfg = [];
        cfg.latency = [-0.8 0.0];
        dataPre = ft_selectdata(cfg, dataCon);
        
        %select stimulation interval
        %cfg = [];
        %cfg.latency = [0.4 1.2];
        %dataPost = ft_selectdata(cfg, selepo);
      
        %% time frequency response with Hann taper with varying length
        %==================================================================
        %mtmconvol gives the time-frequency representation of the trial
        cfg              = [];
        cfg.output       = 'pow';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 7:20;                 
        cfg.toi          = -1:0.01:1.4; 
        cfg.t_ftimwin    = 5./cfg.foi; 
        freq_conv{con}   = ft_freqanalysis(cfg, dataCon); 
        
        % Plot time-frequency for all planar grad 
        cfg = [];
        cfg.parameter    = 'powspctrm';
        cfg.baseline     = [0 inf];
        cfg.baselinetype = 'absolute';
        cfg.showlabels   = 'yes';
        cfg.layout       = 'neuromag306planar.lay';
        
        figure(1); 
        ft_multiplotTFR(cfg, freq_conv{con}); colorbar;
        title (strcat('time-frequency multiplot with Hanning condition', num2str(con)))
        saveas(figure(1),[savepath, subj, '/', subj, '_TF_multiplot_cond_', num2str(con), '.jpeg', ]);
        
        % Plot time-frequency for one channel with max power
        cfg = [];
        cfg.baseline     = [0 inf];
        cfg.baselinetype = 'absolute';
        cfg.maskstyle    = 'saturation';
        cfg.zlim         = 'absmax';
        cfg.channel      = 'MEG2043';
        cfg.layout       = 'neuromag306planar.lay';
        cfg.marker       = 'on';
        
        figure(2);
        subplot(2,1,1); ft_singleplotTFR(cfg, freq_conv{con}); colorbar
        title(strcat('TFR Hanning single channel MEG2043 (grad) for with max power in condition ', num2str(con)));
        
        cfg.layout       = 'neuromag306mag.lay';
        subplot(2,1,2); ft_singleplotTFR(cfg, freq_conv{con}); colorbar
        title(strcat('TFR Hanning single channel MEG2043 (mag) for with max power in condition ', num2str(con)));
        saveas(figure(2),[savepath, subj, '/', subj, '_TF_Hann_single_ch_cond_', num2str(con), '.jpeg', ]);
        
        % Plot topographic representation of the time-frequency representations 
        cfg = [];
        cfg.baseline     = [0.1 inf];
        cfg.baselinetype = 'absolute';
        cfg.xlim         = [-0.8 0];
        cfg.zlim         = 'absmax';
        cfg.ylim         = [7 20];
        cfg.layout       = 'neuromag306planar.lay';
        cfg.marker       = 'on';
        
        figure(3);
        subplot(2,1,1); ft_topoplotTFR(cfg, freq_conv{con}); colorbar
        title (strcat('TFR Hann topoplot (grad) in condition ', num2str(con)))
        
        cfg.layout       = 'neuromag306mag.lay';
        subplot(2,1,2); ft_topoplotTFR(cfg, freq_conv{con}); colorbar
        title(strcat('TFR Hann topoplot (mag) in condition ', num2str(con)));
        saveas(figure(3),[savepath, subj, '/', subj, '_TF_Hann_topo_cond', num2str(con), '.jpeg', ]);
        

        %% TFR with Morlet wavelets 
        %==================================================================
        % To enhance the sensitivity to the non-phase locked activity, we can subtract the phase locked activity 
        
        % Get timelocked 
        timelock = ft_timelockanalysis([], dataCon)
        
        epochs2  = dataCon;
        for i = 1:length(dataCon)
            epochs2.trial{i} = dataCon.trial{i} - timelock.avg;
        end

        % TFR with Morlet wavelets
        cfg = [];
        cfg.method       = 'wavelet';
        cfg.foi          = 7:20;    
        cfg.toi          = -1:0.01:1.4; 
        cfg.width        = 7;
        cfg.pad          = 'nextpow2';
        tfr_wavelet{con} = ft_freqanalysis(cfg, epochs2);
        
        % Plot time-frequency for all planar grad 
        cfg = [];
        cfg.parameter    = 'powspctrm';
        cfg.layout       = 'neuromag306planar';
        cfg.showlabels   = 'yes';
        cfg.baselinetype = 'absolute';
        cfg.baseline     = [0 inf];     
        
        figure(4); ft_multiplotTFR(cfg, tfr_wavelet{con}); colorbar; % wavelet analysis
        title (strcat('TFR with Morlet wavelets in condition', num2str(con)))
        saveas(figure(4),[savepath, subj, '/', subj, '_TFR_wavelets_cond_', num2str(con), '.jpeg']);
        
        % Plot TFR wavelet for one channel with max power
        cfg = [];
        cfg.baseline     = [0 inf];
        cfg.baselinetype = 'absolute';
        cfg.maskstyle    = 'saturation';
        cfg.zlim         = 'absmax';
        cfg.channel      = 'MEG2043';
        cfg.layout       = 'neuromag306planar.lay';
        cfg.marker       = 'on';
        
        figure(5);
        subplot(2,1,1); ft_singleplotTFR(cfg, tfr_wavelet{con}); colorbar
        title(strcat('single channel wavelet MEG2043 (grad) with max power in condition ', num2str(con)));
        
        cfg.layout       = 'neuromag306mag.lay';
        subplot(2,1,2); ft_singleplotTFR(cfg, tfr_wavelet{con}); colorbar
        title(strcat('single channel wavelet MEG2043 (mag) with max power in condition ', num2str(con)));
        saveas(figure(5),[savepath, subj, '/', subj, '_TFwavelet_single_ch_cond_', num2str(con), '.jpeg', ]);
        
        % Plot topography of the time-frequency representations wavelet
        cfg = [];
        cfg.baseline     = [0.1 inf];
        cfg.baselinetype = 'absolute';
        cfg.xlim         = [-0.8 0];
        cfg.zlim         = 'absmax';
        cfg.ylim         = [7 20];
        cfg.layout       = 'neuromag306planar.lay';
        cfg.marker       = 'on';
        
        figure(6);
        subplot(2,1,1); ft_topoplotTFR(cfg, tfr_wavelet{con}); colorbar
        title (strcat('TFR Wavelet topoplot (grad) in condition ', num2str(con)))
        
        cfg.layout       = 'neuromag306mag.lay';
        subplot(2,1,2); ft_topoplotTFR(cfg, tfr_wavelet{con}); colorbar
        title(strcat('TFR Wavelet topoplot (mag) in condition ', num2str(con)));
        saveas(figure(6),[savepath, subj, '/', subj, '_TF_Wavelet_topo_cond', num2str(con), '.jpeg', ]);
        
        %% do spectral analysis using multitapers
        %==================================================================
        % This can be later used to extract individual alpha peak frequency (IAF).
        
        %Get PSD with multitapers
        cfg = [];
        cfg.method       = 'mtmfft';
        cfg.output       = 'pow'; % 'fourier'
        cfg.taper        = 'dpss'; %multitapers
        cfg.keeptrials   = 'yes';
        cfg.pad          = 10; % to [vertually] increase spectral resolution
        cfg.tapsmofrq    = 2;
        cfg.foilim       = [1 50]; %alpharange; %freq band of interest
        freqPre_fft{con} = ft_freqanalysis(cfg, dataPre);  % trials x Ch x freq
        %freqPost        = ft_freqanalysis(cfg, dataPost);
         
        %average over epochs and frequency bins of the alpha range
        cfg=[];
        cfg.frequency    = [7 14]; 
        cfg.avgoverrpt   = 'yes';
        cfg.avgoverfreq  ='yes';
        freqPre_avg{con} = ft_selectdata (cfg, freqPre_fft{con});
        
        freqPre_avg{con}.powspctrm = freqPre_avg{con}.powspctrm'; 
        MAXmag{con} = max(freqPre_avg{con}.powspctrm(1:3:306));  % for plot, to use the same scale for all events
        MAXgrad{con} = max(freqPre_avg{con}.powspctrm(3:3:306));
        
        % Plot topography of the power spectrum
        cfg=[];
        cfg.layout       = 'neuromag306planar.lay'; % plot only gradiometers
        cfg.zlim         = [-1*MAXgrad{con}, MAXgrad{con}];
        cfg.colorbar     = 'yes';
        figure(7); 
        subplot(2,1,1); ft_topoplotER(cfg, freqPre_avg{con});
        title (strcat('PSD for grad with mutitapers in condition', num2str(con)))
        
        cfg.layout       = 'neuromag306mag.lay'; % plot only magnetometers
        cfg.zlim         = [-1*MAXmag{con}, MAXmag{con}];  
        subplot(2,1,2); ft_topoplotER(cfg, freqPre_avg{con});
        title (strcat('PSD for mag with mutitapers in condition', num2str(con)))
        saveas(figure(7),[savepath, subj, '/', subj, '_PSD_multitapers_cond_', num2str(con), '.jpeg']);
        %% do spectral analysis using Hanning taper, for freq < 30 Hz
        %==================================================================
        %Get PSD with Hanning taper
        
        cfg = [];
        cfg.method       = 'mtmfft';
        cfg.output       = 'pow'; 
        cfg.taper        = 'hanning'; %Hanning taper
        cfg.keeptrials   = 'yes';
        cfg.foilim       = [1 30]; %alpharange; 
        freqPre_fft_hann{con} = ft_freqanalysis(cfg, dataPre);  %trials x Ch x freq

        % here we average over epochs and frequency bins of the alpha range
        cfg=[];
        cfg.frequency    = [7 14]; % Hz
        cfg.avgoverrpt   = 'yes';
        cfg.avgoverfreq  ='yes';
        freqPre_avg_hann{con} = ft_selectdata (cfg, freqPre_fft_hann{con});

        freqPre_avg_hann{con}.powspctrm = freqPre_avg_hann{con}.powspctrm'; 

        MAXmag_hann{con}  = max(freqPre_avg_hann{con}.powspctrm(1:3:306));  % for plot, to use the same scale for all events
        MAXgrad_hann{con} = max(freqPre_avg_hann{con}.powspctrm(3:3:306));

        % Plot topography of the power spectrum
        cfg=[];
        cfg.layout = 'neuromag306planar.lay';
        cfg.zlim = [-1*MAXgrad_hann{con}, MAXgrad_hann{con}];
        cfg.colorbar = 'yes';
        figure(7); 
        subplot(2,1,1); ft_topoplotER(cfg, freqPre_avg_hann{con});
        title (strcat('PSD for grad with Hanning taper in condition', num2str(con)));
        
        cfg.layout = 'neuromag306mag.lay'; 
        cfg.zlim = [-1*MAXmag_hann{con}, MAXmag_hann{con}];
        cfg.colorbar = 'yes';
        subplot(2,1,2); ft_topoplotER(cfg, freqPre_avg_hann{con});
        title (strcat('PSD for mag with Hanning taper in condition', num2str(con)));
        saveas(figure(7),[savepath, subj, '/', subj, '_PSD_hann_cond_', num2str(con), '.jpeg']);
        
    end
    
    %% Compute contrast between two conditions for TFR with Hanning taper
    %======================================================================
    cfg            = [];
    cfg.parameter  = 'powspctrm';
    cfg.operation  = '(x1-x2)/(x1+x2)';
    diff_hann_v3v1 = ft_math(cfg, freq_conv{3}, freq_conv{1});
    
    % Plot topography of the difference
    cfg            = [];
    cfg.layout     = 'neuromag306planar.lay';
    cfg.marker     = 'on';
    figure(8);
    subplot(2,1,1); ft_topoplotTFR(cfg, diff_hann_v3v1); colorbar;
    title ('Diff topo (grad) after TFR with Hanning between fast and slow condition')
    
    cfg.layout     = 'neuromag306mag.lay';
    subplot(2,1,2); ft_topoplotTFR(cfg, diff_hann_v3v1); colorbar;
    title ('Diff topo (mag) after TFR with Hanning between fast and slow condition')
    saveas(figure(8),[savepath, subj, '/', subj, '_Diff_v3v1_TFR_Hann.jpeg']);
    
    %% Compute contrast between two conditions for TFR with Morlet wavelets 
    %======================================================================
    cfg            = [];
    cfg.parameter  = 'powspctrm';
    cfg.operation  = '(x1-x2)/(x1+x2)';
    diff_wavelet_v3v1 = ft_math(cfg, tfr_wavelet{3}, tfr_wavelet{1});
    
    % Plot topography of the difference
    cfg            = [];
    cfg.layout     = 'neuromag306planar.lay';
    cfg.marker     = 'on';
    figure(9);
    subplot(2,1,1); ft_topoplotTFR(cfg, diff_hann_v3v1); colorbar;
    title ('Diff topo (grad) after TFR with wavelets between fast and slow condition')
    
    cfg.layout     = 'neuromag306mag.lay';
    subplot(2,1,2); ft_topoplotTFR(cfg, diff_wavelet_v3v1); colorbar;
    title ('Diff topo (mag) after TFR with wavelets between fast and slow condition')
    saveas(figure(9),[savepath, subj, '/', subj, '_Diff_v3v1_TFR_wavelet.jpeg']);
    
    %% Calculate absolute power differences between conditions for PSD with multitapers
    %======================================================================
    
    diff_v3v1 = freqPre_avg{1}; % use freqPre{1} structure, but change power for the 'difference'
    diff_v3v1.powspctrm = freqPre_avg{3}.powspctrm  - freqPre_avg{1}.powspctrm;
    
    % Plot topography of the difference
    cfg=[];
    cfg.layout     = 'neuromag306planar.lay'; % plot only grad
    cfg.colorbar   = 'yes';
    figure(10); 
    subplot(2,1,1); ft_topoplotER(cfg, diff_v3v1);
    title ('Diff PSD with multitapers (grad) Fast-Slow');
    
    cfg.layout     = 'neuromag306mag.lay'; % plot only mag
    cfg.colorbar   = 'yes';
    subplot(2,1,2); ft_topoplotER(cfg, diff_v3v1);
    title ('Diff PSD with multitapers (mag) Fast-Slow');
    saveas(figure(10),[savepath, subj, '/', subj, 'Diff_v3v1_PSD_multitapers.jpeg']);

    %% Calculate absolute power differences between conditions for PSD with Hanning taper
    %======================================================================
    
    diff_v3v1_hann = freqPre_avg_hann{1}; % use freqPre{1} structure, but change power for the 'difference'
    diff_v3v1_hann.powspctrm = freqPre_avg_hann{3}.powspctrm  - freqPre_avg_hann{1}.powspctrm;
    
    % Plot topography of the difference
    cfg=[];
    cfg.layout = 'neuromag306planar.lay'; % plot only gradiometers
    cfg.colorbar = 'yes';
    figure(11); 
    subplot(2,1,1); ft_topoplotER(cfg, diff_v3v1_hann);
    title ('Diff PSD with Hanning taper (grad) Fast-Slow');
    
    cfg.layout = 'neuromag306mag.lay'; % plot only magnetometers
    cfg.colorbar = 'yes';
    subplot(2,1,2); ft_topoplotER(cfg, diff_v3v1_hann);
    title ('Diff PSD with Hanning taper (mag) Fast-Slow');
    saveas(figure(11),[savepath, subj, '/', subj, 'Diff_v3v1_PSD_hann_.jpeg']);
    
    %% Averaged time-frequency responses. hanning with fixed window.
    
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmconvol';
    cfg.taper        = 'hanning';
    cfg.toi          = -0.8 : 0.1 : 0.0;
    cfg.foi          = 7:20;
    cfg.t_ftimwin    = ones(size(cfg.foi)) * 0.5;
    TFR_all          = ft_freqanalysis(cfg, epochs);

    cfg.trials       = find(allinfo.prev_stim_type==2);
    TFR_slow         = ft_freqanalysis(cfg, epochs);

    cfg.trials       = find(allinfo.prev_stim_type==8);
    TFR_fast         = ft_freqanalysis(cfg, epochs);
 
    cfg = [];
    cfg.baseline     = [0 inf];
    cfg.baselinetype = 'absolute';
    cfg.layout       = 'neuromag306mag.lay';

    figure(1);
    ft_multiplotTFR(cfg, TFR_slow); 
    title('TFR Hanning with fxed window in slow condition'); colorbar;
    saveas(figure(1),[savepath, subj, '/', subj, '_TFR_Hann_fixwindow_slow.jpeg' ]);
    
    figure(2);
    ft_multiplotTFR(cfg, TFR_fast);
    title('TFR Hanning with fxed window in fast condition'); colorbar;
    saveas(figure(2),[savepath, subj, '/', subj, '_TFR_Hann_fixwindow_fast.jpeg' ]);
    
    %Compute contrast between response for different velocity 
    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = '(x1-x2)/(x1+x2)';
    TFR_diff      = ft_math(cfg, TFR_fast, TFR_slow);

    cfg = [];
    cfg.marker  = 'on';
    cfg.layout  = 'neuromag306mag.lay';
    figure; ft_multiplotTFR(cfg, TFR_diff);

end