% ft_pms_step1_sensors.m
%
% Load data epochs 
% Load info about the preseeding trials
% Calculate spectral power on sensors in [-.8 to 0] interval
% plot average alpha power  distribution of epochs following different
% types of stimuli (Slow, Medium, Fast); plot Fast-Slow power differences
% save pictures as PPTX

%restoredefaultpath;
clear;
close all;
clc;
screensize = get( groot, 'Screensize' );

tapsmofrq  = 2;
megfolder = strcat( 'meg_sensors_tapsmofrq', num2str(tapsmofrq));
alpharange = [7 14];

%%
% NB: add to path: 
fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(path, fieldtripfolder)
ft_defaults;
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
    
    %% Plot average ERF: this has no sense for rebound, but you will see the evoked response
    cfg = [];
    cfg.channel = epochs.label;
    cfg.covariance = 'yes';
    avg = ft_timelockanalysis(cfg,epochs);
    
    hh=figure;
    ttt = find(avg.time>-0.3 & avg.time<0.5);
    plot (avg.time(ttt), avg.avg(:, ttt));
    title ('Sensor averages')
   
    %% Calculating power in sensors, find max weighted freq, freq range and P of gamma enhancement in this freq range.
    hh=figure;
    LineColor{1} = 'b'; LineColor{2} = 'g'; LineColor{3} = 'r';
    LineColorMax{1} = '--b'; LineColorMax{2} = '--g'; LineColorMax{3} = '--r';
    ax=gca;
    
%     % posterior sensors, you may use them later for extracting alpha peak
%     power
%     Ch = {'MEG1932',  'MEG1922', 'MEG2042',  'MEG2032',  'MEG2112', 'MEG2122',  'MEG2342', 'MEG2332',  'MEG1732', 'MEG1942', 'MEG1912', 'MEG2012', 'MEG2022', 'MEG2312', 'MEG2322', 'MEG2512',...
%           'MEG1933',  'MEG1923', 'MEG2043',  'MEG2033',  'MEG2113', 'MEG2123',  'MEG2343', 'MEG2333',  'MEG1733', 'MEG1943', 'MEG1913', 'MEG2013', 'MEG2023', 'MEG2313', 'MEG2323', 'MEG2513'};
   
    %%    
    for con=1:3 % for conditions
        %select epochs according to preceding trials
        cfg = [];
        cfg.trials = prev{con}; % EV{con}';
        selepo = ft_selectdata(cfg, epochs);
        
        cfg = [];
        cfg.latency = [-0.8 0.0];
        dataPre = ft_selectdata(cfg, selepo);
        
        
% %         cfg.latency = [0.4 1.2];
% %         latency_post = cfg.latency;
% %         [dataPost] = ft_selectdata(cfg, selepo);
        
        % do spectral analysis using multitapers, we use 'pad' to make smooth
        % spectrum. This can be later used to extract individual alpha peak frequency (IAF).
        cfg = [];
        cfg.method = 'mtmfft';
        cfg.output ='pow'; % 'fourier'
        cfg.taper = 'hanning';
        cfg.keeptrials = 'yes';
        cfg.foilim = alpharange; %freq band of interest
        freqPre{con} = ft_freqanalysis(cfg, selepo);  % trials x Ch x freq
        %freqPost = ft_freqanalysis(cfg, dataPost);
        
        cfg              = [];
        cfg.output       = 'pow';
        cfg.channel      = 'all';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 7:0.5:15;                         % analysis 2 to 30 Hz in steps of 2 Hz
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
        cfg.toi          = -1:0.05:1.2; 
        
        cfg = [];
        cfg.baseline     = [0.2 1.2];
        cfg.baselinetype = 'absolute';
        cfg.showlabels   = 'yes';
        cfg.layout       = 'neuromag306planar.lay';
        figure; ft_multiplotTFR(cfg, freqPre{con});

        
        cfg = [];
        cfg.baseline     = [0.2 1.2];
        cfg.baselinetype = 'absolute';
        cfg.maskstyle    = 'saturation';
        %cfg.zlim         = [0 2.33e-23];
        cfg.channel      = 'MEG2042';
        cfg.layout       = 'neuromag306planar.lay';
        
        figure;
        ft_singleplotTFR(cfg, freqPre{con}); title('fast');

        % here we average over epochs and frequency bins of the alpha range
        cfg=[];
        cfg.frequency = alpharange; % Hz
        cfg.avgoverrpt = 'yes';
        %cfg.avgoverfreq ='yes';
        freqPre{con} = ft_selectdata (cfg, freqPre{con});
        
        freqPre{con}.powspctrm = freqPre{con}.powspctrm';
        
        %MAXmag(con) = max(freqPre.powspctrm(1:3:306))  % for plot, to use the same scale for all events
        MAXgrad(con) = max(freqPre{con}.powspctrm(3:3:204));  % for plot, to use the same scale for all events

    end
    %% calculate absolute power differences  between conditions: afterFast - afterSlow
    diff_v3v1=freqPre{1} ; % use freqPre{1} structure, but change power for the 'difference'
    diff_v3v1.powspctrm = freqPre{3}.powspctrm  - freqPre{1}.powspctrm;

    %% plot power and power difference (Fast-Slow) for gradiometers
    %pos_fig1 = [10 10 screensize(3)/10*9 screensize(4)/3];    
    
    pos_fig1 = [10 10 screensize(3)/10*9 screensize(4)/3];    
    hh=figure('Position',pos_fig1);

    cfg=[];
    cfg.layout = 'neuromag306planar.lay'; % plot only gradeometers
    cfg.zlim = [-1*max(MAXgrad), max(MAXgrad)];
    cfg.colorbar = 'yes';
    subplot(1,4,1)
    ft_topoplotER(cfg,freqPre{1}); title('following Slow');
    subplot(1,4,2)
    ft_topoplotER(cfg,freqPre{2}); title('following Medium');
    subplot(1,4,3)
    ft_topoplotER(cfg,freqPre{3}); title('following Fast');
    cfg=[];
    cfg.layout = 'neuromag306planar.lay'; % plot only gradeometers
    subplot(1,4,4); ft_topoplotER(cfg,diff_v3v1); title('alpha power Fast-Slow');
    
    subtitle ('Gradiometers Alpha power [-0.8 to 0 s], according to the type of the preceeding event')
    

%     %% plot power  and power difference (Fast-Slow) for magnetometers
%     pos_fig1 = [10 10 screensize(3)/10*9 screensize(4)/3];    
%     hh=figure('Position',pos_fig1);
%     cfg=[];
%     cfg.layout = 'neuromag306mag.lay'; % plot only magnetometers
%     cfg.zlim = [-1*max(MAXmag), max(MAXmag)];
%     subplot(1,4,1); ft_topoplotER(cfg,freqPre{1}); title('following Slow');
%     subplot(1,4,2); ft_topoplotER(cfg,freqPre{2}); title('following Medium');
%     subplot(1,4,3); ft_topoplotER(cfg,freqPre{3}); title('following Fast');
%     cfg=[];
%     cfg.layout = 'neuromag306mag.lay'; % plot only magnetometers
%     subplot(1,4,4); ft_topoplotER(cfg,diff_v3v1); title('alpha power Fast-Slow');
% %     colorbar
% %     
    

end % for subjects


