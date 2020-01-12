%IRASA
% the extraction of rhythmic spectral features 
% from the electrophysiological signal based on Irregular Resampling
% Auto-Spectral Analysis (IRASA, Wen & Liu, Brain Topogr. 2016)

clear all;
close all;

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(path, fieldtripfolder)
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);

savepath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';
realdatapath = '/home/a_shishkina/data/KI/SUBJECTS/';

%load subj list
     
SUBJ_NT = [ '0101'; '0102'; '0103'; '0104'; '0105'; '0135'; '0136';...  
            '0137'; '0138'; '0140'; '0158'; '0162'; '0163'; '0178';...
            '0179'; '0255'; '0257'; '0348'; '0378'; '0379'; '0384']; 
        
SUBJ_ASD = ['0106'; '0107'; '0139'; '0141'; '0159'; '0160'; '0161';...  
            '0164'; '0253'; '0254'; '0256'; '0273'; '0274'; '0275';...
            '0276'; '0346'; '0347'; '0349'; '0351'; '0358';...
            '0380'; '0381'; '0382'; '0383'];  
%without '0357';
SUBJ = [SUBJ_ASD; SUBJ_NT];

%load posterior sensors list
post_sens = {'MEG1932',  'MEG1922', 'MEG2042',  'MEG2032',  'MEG2112', 'MEG2122',  'MEG2342', 'MEG2332',  'MEG1732', 'MEG1942', 'MEG1912', 'MEG2012', 'MEG2022', 'MEG2312', 'MEG2322', 'MEG2512',...
             'MEG1933',  'MEG1923', 'MEG2043',  'MEG2033',  'MEG2113', 'MEG2123',  'MEG2343', 'MEG2333',  'MEG1733', 'MEG1943', 'MEG1913', 'MEG2013', 'MEG2023', 'MEG2313', 'MEG2323', 'MEG2513'};

SUBJ = ['0384'];
%loop for all subjects
for s=1: size (SUBJ,1)
    
    subj = SUBJ (s,:); 
    savemegto = strcat(savepath, subj);
    epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
    
    epo = load(strcat(epofolder, subj, '_preproc_epochs.mat'));
    
    %select grad and mag epochs separately for slow and fast conditions
    cfg = [];
    cfg.channel = 'MEGMAG';
    cfg.latency = [-0.8 0.0];
    data_fast_grad = ft_selectdata(cfg, epo.fast_epochs);
    data_slow_grad = ft_selectdata(cfg, epo.slow_epochs);
    
    cfg.channel = 'MEGGRAD';
    cfg.latency = [-0.8 0.0];
    data_fast_mag = ft_selectdata(cfg, epo.fast_epochs);
    data_slow_mag = ft_selectdata(cfg, epo.slow_epochs);
    
    % perform IRASA and regular spectral analysis
    epochs = [data_slow_grad; data_fast_grad; data_slow_mag; data_fast_mag];
    for s = 1: size(epochs,1)
        data = epochs(s,:);
        
        cfg               = [];
        cfg.method        = 'irasa';
        cfg.taper         = 'hanning'; 
        cfg.pad           = 1024/1000; %sampl/fs
        cfg.foilim        = [2 40];       
        cfg.tapsmofrq     = 3; 
        cfg.keeptrials    = 'yes';
        frac_r{s} = ft_freqanalysis(cfg, data);
        
        cfg.method     = 'mtmfft';
        orig_r{s} = ft_freqanalysis(cfg, data);
        
        cfg = [];
        cfg.avgoverchan = 'yes';
        frac_avg{s} = ft_selectdata(cfg, frac_r{s});
        orig_avg{s} = ft_selectdata(cfg, orig_r{s});
        
        % average across the sub-segments
        frac_s = {}; 
        orig_s = {};
        for rpt = unique(frac_avg{s}.trialinfo)'
            cfg               = [];
            cfg.trials        = find(frac_avg{s}.trialinfo==rpt);
            cfg.avgoverrpt    = 'yes';
            frac_s{end+1} = ft_selectdata(cfg, frac_avg{s});
            orig_s{end+1} = ft_selectdata(cfg, orig_avg{s});
        end
        frac_a{s} = ft_appendfreq([], frac_s{:});
        orig_a{s} = ft_appendfreq([], orig_s{:});

        % average across trials
        cfg               = [];
        cfg.trials        = 'all';
        cfg.avgoverrpt    = 'yes';
        frac{s} = ft_selectdata(cfg, frac_a{s});
        orig{s} = ft_selectdata(cfg, orig_a{s});

        % subtract the fractal component from the power spectrum
        cfg               = [];
        cfg.parameter     = 'powspctrm';
        cfg.operation     = 'x2-x1';
        osci{s} = ft_math(cfg, frac{s}, orig{s});
        
        
        % plot the fractal component and the power spectrum 
        if s == 1 | s == 2
            figure(1);
            subplot(2,1,s);
            plot(frac{s}.freq, frac{s}.powspctrm, ...
              'linewidth', 2, 'color', [0 0 0])
            hold on;
            plot(orig{s}.freq, orig{s}.powspctrm, ...
              'linewidth', 2, 'color', [0, 0.4470, 0.7410])
            xlim([2,40]);
            hold off; 

            title_name = ['after slow, grad'; 'after fast, grad'; 'after slow,  mag'; 'after fast,  mag'];
          
            % plot the full-width half-maximum of the oscillatory component
            f = fit(osci{s}.freq', osci{s}.powspctrm(1,:)', 'gauss1');
            mean = f.b1;
            std  = f.c1/sqrt(2)*2.3548;
            fwhm = [mean-std/2 mean+std/2];
            yl   = get(gca, 'YLim');
            p = patch([fwhm flip(fwhm)], [yl(1) yl(1) yl(2) yl(2)], [1 1 1]);
            uistack(p, 'bottom');  title([subj, ', IRASA power, ' title_name(s,:)])
            legend('FWHM oscillation', 'Fractal component', 'Power spectrum');
            xlabel('Frequency'); ylabel('Power');
            set(gca, 'YLim', yl); 
        else
            figure(2);
            subplot(2,1,s-2);
            plot(frac{s}.freq, frac{s}.powspctrm, ...
              'linewidth', 2, 'color', [0 0 0])
            hold on; plot(orig{s}.freq, orig{s}.powspctrm, ...
              'linewidth', 2, 'color', [0.8500, 0.3250, 0.0980])
            xlim([2,40]);
            hold off; 

            % plot the full-width half-maximum of the oscillatory component
            f    = fit(osci{s}.freq', osci{s}.powspctrm(1,:)', 'gauss1');
            mean = f.b1;
            std  = f.c1/sqrt(2)*2.3548;
            fwhm = [mean-std/2 mean+std/2];
            yl   = get(gca, 'YLim');
            p = patch([fwhm flip(fwhm)], [yl(1) yl(1) yl(2) yl(2)], [1 1 1]);
            uistack(p, 'bottom');  title([subj, ', IRASA power, ' title_name(s,:)])
            legend('FWHM oscillation', 'Fractal component', 'Power spectrum');
            xlabel('Frequency'); ylabel('Power');
            set(gca, 'YLim', yl); 
        end    
    end
    %save figures
    saveas(figure(1), [savepath, subj, '/', subj, '_IRASA_power_grad.jpeg']);
    saveas(figure(2), [savepath, subj, '/', subj, '_IRASA_power_mag.jpeg']);
end