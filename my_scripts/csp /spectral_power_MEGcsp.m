tic
%double check to estimate the spectral power in the MEG_csp data for all eigenvalue
for s = 1:6
        for j = 1:size(Xcsp_fast,1)
            trial_fast(j,:,:) = squeeze(Xcsp_fast(j,:,:))*A{s};
            csp_data_fast{j} = squeeze(trial_fast(j,:,:))';
        end
        fast{s} = csp_data_fast;
        
        for j = 1:size(Xcsp_slow,1)     
            trial_slow(j,:,:) = squeeze(Xcsp_slow(j,:,:))*A{s};
            csp_data_slow{j} = squeeze(trial_slow(j,:,:))';
        end
        slow{s} = csp_data_slow;
        
        data_comp_fast{s} = data_fast;
        data_comp_fast{s}.trial = fast{s};

        data_comp_slow{s} = data_slow;
        data_comp_slow{s}.trial = slow{s};
        
        cfg = [];
        cfg.method       = 'mtmfft';
        cfg.output       = 'pow'; 
        cfg.taper        = 'hanning'; 
        cfg.pad          = 10; 
        cfg.foilim       = [5 30];
        cfg.tapsmofrq    = 3; 
        fft_fast  = ft_freqanalysis(cfg, data_comp_fast{s}); 
        fft_slow   = ft_freqanalysis(cfg, data_comp_slow{s});

        cfg = [];
        cfg.avgoverchan = 'yes';
        fft_fast = ft_selectdata(cfg,fft_fast);
        fft_slow = ft_selectdata(cfg,fft_slow);


        figure(1)
        subplot(2,3,s)
        plot(fft_fast.freq, fft_fast.powspctrm, '-r')
        hold on;
        plot(fft_slow.freq, fft_slow.powspctrm, '-b')
        legend('fast', 'slow'); 
        title(['compoment', num2str(s)])
end
toc

%save figures
saveas(figure(1), [savepath, subj, '/', subj, '_spectral_power_MEG_csp.jpeg']);