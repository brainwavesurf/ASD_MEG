% Load data epochs 
% Load info about the preceding trials
% Calculate spectral power on sensors in [-.8 to 0] interval
% Plot the time-frequency representation of the trial for grad in the each condition by mtmconvol 
% Plot the time-frequency representation for one channel with max power according to previous plot
% Plot the time-frequency topoplot for grad in specific time and frequency range  
% Plot TFR (time-frequency response) by Morlet wavelets with subtraction the phase locked activity to enhance the sensitivity to the induced activity 
% Plot PSD (power spectra density) with multitapers
% Do averaging over epochs and frequency bins of the alpha range
% Plot Fast-Slow power differences for grad and mag
% Plot PSD (power spectra density) with Hanning taper for low freq
% Do averaging over epochs and frequency bins of the alpha range
% Plot Fast-Slow power differences for grad and mag
%%
clear all;
close all;

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path('/home/a_shishkina/fieldtrip/external/mne/', path);
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
        selepo = ft_selectdata(cfg, epochs);
        
        %select interstimulus interval
        cfg = [];
        cfg.latency = [-0.8 0.0];
        dataPre = ft_selectdata(cfg, selepo);
        
        %select stimulation interval
        %cfg = [];
        %cfg.latency = [0.4 1.2];
        %dataPost = ft_selectdata(cfg, selepo);
      
        %% mtmconvol gives the time-frequency representation of the trial, how the frequency content changes over time.
        % time frequency response with Hann taper with varying length
        cfg              = [];
        cfg.output       = 'pow';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 7:1:20;                 
        cfg.toi          = -1:0.01:1.4; 
        cfg.t_ftimwin    = 5./cfg.foi; 
        freqPre_conv{con}= ft_freqanalysis(cfg, selepo); 
        
        
        % Plot time-frequency for all planar grad 
        cfg = [];
        cfg.parameter    = 'powspctrm';
        cfg.baseline     = [0.1 inf];
        cfg.baselinetype = 'absolute';
        cfg.showlabels   = 'yes';
        cfg.layout       = 'neuromag306planar.lay';
        %cfg.channel     = post_sens;
        %cfg.masknans     = 'no'
        %cfg.renderer    = 'painters';
        
        figure(1); 
        ft_multiplotTFR(cfg, freqPre_conv{con}); colorbar;
        title (strcat('time-frequency multiplot with Hanning condition', num2str(con)))
        saveas(figure(1),[savepath, subj, '/', subj, '_TF_multiplot_cond_', num2str(con), '.jpeg', ]);
        
        %% Plot time-frequency for one channel with max power
        cfg = [];
        cfg.baseline     = [0.1 inf];
        cfg.baselinetype = 'absolute';
        cfg.maskstyle    = 'saturation';
        cfg.zlim         = [-1.99e-24 2.4e-23];
        cfg.channel      = 'MEG2043';
        cfg.layout       = 'neuromag306planar.lay';
        cfg.marker       = 'on';
        
        figure(2);
        subplot(2,1,1); ft_singleplotTFR(cfg, freqPre_conv{con}); colorbar
        title(strcat('TFR Hanning single channel MEG2043 (grad) for with max power in condition ', num2str(con)));
        
        cfg.layout       = 'neuromag306mag.lay';
        subplot(2,1,2); ft_singleplotTFR(cfg, freqPre_conv{con}); colorbar
        title(strcat('TFR Hanning single channel MEG2043 (mag) for with max power in condition ', num2str(con)));
        saveas(figure(2),[savepath, subj, '/', subj, '_TF_Hann_single_ch_cond_', num2str(con), '.jpeg', ]);
        
        %% A topographic representation of the time-frequency representations 
        cfg = [];
        cfg.baseline     = [0.1 inf];
        cfg.baselinetype = 'absolute';
        cfg.xlim         = [-0.3 0];
        cfg.zlim         = [-5.01e-24 5.01e-24];
        cfg.ylim         = [9.5 13];
        cfg.layout       = 'neuromag306planar.lay';
        cfg.marker       = 'on';
        
        figure(3);
        subplot(2,1,1); ft_topoplotTFR(cfg, freqPre_conv{con}); colorbar
        title (strcat('TFR Hann topoplot (grad) in condition ', num2str(con)))
        
        cfg.zlim         = [-5.01e-27 5.01e-27];
        cfg.layout       = 'neuromag306mag.lay';
        subplot(2,1,2); ft_topoplotTFR(cfg, freqPre_conv{con}); colorbar
        title(strcat('TFR Hann topoplot (mag) in condition ', num2str(con)));
        saveas(figure(3),[savepath, subj, '/', subj, '_TF_Hann_topo_cond', num2str(con), '.jpeg', ]);
        
        
        %% identify the indices of trials with high and low alpha power
        %tmp     = mean(freqPre_conv{con}.powspctrm,2); %mean over frequencies between 7-14Hz
        %ind     = find(mean(tmp,1)==max(mean(tmp,1)));  % find the sensor where power is max
        %indhigh = find(tmp(:,ind)>=median(tmp(:,ind)));   
        
        %% TFR with Morlet wavelets 
        % To enhance the sensitivity to the non-phase locked activity, we can subtract the phase locked activity 
        
        % Get timelocked again
        timelock = ft_timelockanalysis([], selepo)
        epochs2 = selepo;
        for i = 1:length(selepo)
            epochs2.trial{i} = selepo.trial{i} - timelock.avg;
        end

        % TFR with Morlet wavelets
        cfg = [];
        cfg.method       = 'wavelet';
        cfg.foi          = 7:20;    
        cfg.toi          = -1:0.01:1.4; 
        cfg.width        = 7;
        cfg.pad          = 'nextpow2';
        tfr_wavelet{con} = ft_freqanalysis(cfg, epochs2);
        
        cfg = [];
        cfg.parameter    = 'powspctrm';
        cfg.layout       = 'neuromag306planar';
        cfg.showlabels   = 'yes';
        cfg.baselinetype = 'absolute';
        cfg.baseline     = [0.1 inf];     
        
        figure(4); ft_multiplotTFR(cfg, tfr_wavelet{con}); colorbar; % wavelet analysis
        title (strcat('TFR with Morlet wavelets in condition', num2str(con)))
        saveas(figure(4),[savepath, subj, '/', subj, '_TFR_wavelets_cond_', num2str(con), '.jpeg']);
        
        cfg = [];
        cfg.baseline     = [0.1 inf];
        cfg.baselinetype = 'absolute';
        cfg.maskstyle    = 'saturation';
        cfg.zlim         = [-7.02e-22 5.6e-21];
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
        
        %% do spectral analysis using multitapers, we use 'pad' to make smooth
        % spectrum. This can be later used to extract individual alpha peak frequency (IAF).
        
        %Get PSD with multitapers
        cfg = [];
        cfg.method       = 'mtmfft';
        cfg.output       ='pow'; % 'fourier'
        cfg.taper        = 'dpss'; %multitapers
        cfg.keeptrials   = 'yes';
        cfg.pad          = 10; % to [vertually] increase spectral resolution
        cfg.tapsmofrq    = tapsmofrq;
        cfg.foilim       = [1 50]; %alpharange; %freq band of interest
        freqPre_fft{con} = ft_freqanalysis(cfg, dataPre);  % trials x Ch x freq
        %freqPost        = ft_freqanalysis(cfg, dataPost);
         
        % here we average over epochs and frequency bins of the alpha range
        cfg=[];
        cfg.frequency    = [7 14]; % Hz
        cfg.avgoverrpt   = 'yes';
        cfg.avgoverfreq  ='yes';
        freqPre_avg{con} = ft_selectdata (cfg, freqPre_fft{con});
        
        freqPre_avg{con}.powspctrm = freqPre_avg{con}.powspctrm'; 
        
        MAXmag{con} = max(freqPre_avg{con}.powspctrm(1:3:306));  % for plot, to use the same scale for all events
        MAXgrad{con} = max(freqPre_avg{con}.powspctrm(3:3:306));
      
        cfg=[];
        cfg.layout = 'neuromag306planar.lay'; % plot only gradeometers
        cfg.zlim = [-1*MAXgrad{con}, MAXgrad{con}];
        cfg.colorbar = 'yes';
        
        figure(5); 
        subplot(2,1,1); ft_topoplotER(cfg, freqPre_avg{con});
        title (strcat('PSD for grad with mutitapers in condition', num2str(con)))
        
        cfg=[];
        cfg.layout = 'neuromag306mag.lay'; % plot only magnetometers
        cfg.zlim = [-1*MAXmag{con}, MAXmag{con}];
        cfg.colorbar = 'yes';
        
        subplot(2,1,2); ft_topoplotER(cfg, freqPre_avg{con});
        title (strcat('PSD for mag with mutitapers in condition', num2str(con)))
        saveas(figure(5),[savepath, subj, '/', subj, '_PSD_cond_', num2str(con), '.jpeg']);

        
    end
    
    %% diff for mtmconvol
    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = '(x1-x2)/(x1+x2)';

    TFR_diff_v3v1 = ft_math(cfg, tfr_wavelet{3}, tfr_wavelet{1});
    
    cfg = [];
    %cfg.baseline     = [0.1 inf];
    %cfg.baselinetype = 'absolute';
    cfg.xlim         = [-0.3 0];
    %cfg.zlim         = [-1.99e-24 2.4e-23];
    cfg.ylim         = [9.5 11.5];
    cfg.layout       = 'neuromag306planar.lay';
    cfg.marker       = 'on';
    %cfg.channel      = 'MEG*1';

    figure(10);
    ft_topoplotTFR(cfg, TFR_diff_v3v1); colorbar
    title ('time-frequency topoplot after wavelet for diff between v3 and v1 ')
    saveas(figure(10),[savepath, subj, '/', subj, '_TFwavelet_topoplot_diff_v3_v1.jpeg']);
    
    %% calculate absolute power differences  between conditions: afterFast - afterSlow
    diff_v3v1=freqPre_avg{1}; % use freqPre{1} structure, but change power for the 'difference'
    diff_v3v1.powspctrm = freqPre_avg{3}.powspctrm  - freqPre_avg{1}.powspctrm;

    %% plot power and power difference (Fast-Slow) for gradiometers    
    cfg=[];
    cfg.layout = 'neuromag306planar.lay'; % plot only gradiometers
    cfg.colorbar = 'yes';
    
    figure(6); 
    subplot(2,1,1); ft_topoplotER(cfg, diff_v3v1);
    title ('PSD for grad alpha power Fast-Slow');
    
    cfg=[];
    cfg.layout = 'neuromag306mag.lay'; % plot only magnetometers
    cfg.colorbar = 'yes';

    subplot(2,1,2); ft_topoplotER(cfg, diff_v3v1);
    title ('PSD for mag alpha power Fast-Slow');
    saveas(figure(6),[savepath, subj, '/', subj, '_PSD_diff_.jpeg']);

        
    for con=1:3
        cfg = [];
        cfg.trials = prev{con}; % EV{con}';
        selepo = ft_selectdata(cfg, epochs);

        %select interstimulus interval
        cfg = [];
        cfg.latency = [-0.8 0.0];
        dataPre = ft_selectdata(cfg, selepo);
        %% do spectral analysis using Hanning taper, for low freq

        %Get PSD with Hanning taper
        cfg = [];
        cfg.method       = 'mtmfft';
        cfg.output       ='pow'; 
        cfg.taper        = 'hanning'; %Hanning taper
        cfg.keeptrials   = 'yes';
        cfg.foilim       = [1 30]; %alpharange; %freq band of interest
        freqPre_fft_h{con} = ft_freqanalysis(cfg, dataPre);  % trials x Ch x freq
        %freqPost        = ft_freqanalysis(cfg, dataPost);

        % here we average over epochs and frequency bins of the alpha range
        cfg=[];
        cfg.frequency    = [7 14]; % Hz
        cfg.avgoverrpt   = 'yes';
        cfg.avgoverfreq  ='yes';
        freqPre_avg{con} = ft_selectdata (cfg, freqPre_fft_h{con});

        freqPre_avg{con}.powspctrm = freqPre_avg{con}.powspctrm'; 

        MAXmag{con} = max(freqPre_avg{con}.powspctrm(1:3:306));  % for plot, to use the same scale for all events
        MAXgrad{con} = max(freqPre_avg{con}.powspctrm(3:3:306));

        cfg=[];
        cfg.layout = 'neuromag306planar.lay'; % plot only gradeometers
        cfg.zlim = [-1*MAXgrad{con}, MAXgrad{con}];
        cfg.colorbar = 'yes';

        figure(7); 
        subplot(2,1,1); ft_topoplotER(cfg, freqPre_avg{con});
        title (strcat('PSD for grad with Hanning taper in condition', num2str(con)))

        cfg=[];
        cfg.layout = 'neuromag306mag.lay'; % plot only magnetometers
        cfg.zlim = [-1*MAXmag{con}, MAXmag{con}];
        cfg.colorbar = 'yes';

        subplot(2,1,2); ft_topoplotER(cfg, freqPre_avg{con});
        title (strcat('PSD for mag with single Hanning taper in condition', num2str(con)))
        saveas(figure(7),[savepath, subj, '/', subj, '_PSD_hann_cond_', num2str(con), '.jpeg']);
    end
    %% calculate absolute power differences  between conditions: afterFast - afterSlow
    diff_v3v1=freqPre_avg{1}; % use freqPre{1} structure, but change power for the 'difference'
    diff_v3v1.powspctrm = freqPre_avg{3}.powspctrm  - freqPre_avg{1}.powspctrm;

    %% plot power and power difference (Fast-Slow) for gradiometers    
    cfg=[];
    cfg.layout = 'neuromag306planar.lay'; % plot only gradiometers
    cfg.colorbar = 'yes';
    
    figure(8); 
    subplot(2,1,1); ft_topoplotER(cfg, diff_v3v1);
    title ('PSD Hann for grad alpha power Fast-Slow');
    
    cfg=[];
    cfg.layout = 'neuromag306mag.lay'; % plot only magnetometers
    cfg.colorbar = 'yes';

    subplot(2,1,2); ft_topoplotER(cfg, diff_v3v1);
    title ('PSD Hann for mag alpha power Fast-Slow');
    saveas(figure(8),[savepath, subj, '/', subj, '_PSD_hann_diff_.jpeg']);
end % for subjects