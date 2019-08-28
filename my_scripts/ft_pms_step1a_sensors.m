% ft_pms_step1_sensors.m
% Do TMF on 4 max sensors and save sensor-level parameters in subject's folder
% (wF, p, etc.)

%restoredefaultpath;
clear;
close all;
clc;

tapsmofrq  = 5;
gridres = 6; % 6 mm grid
megfolder = strcat( 'meg_sensors_tapsmofrq', num2str(tapsmofrq));
gammarange = [35 110];
lambda=5;

%%
% NB: add to path: 
fieldtripfolder = '/home/kolai/Documents/Shishkina/ProgramFiles/Fieldtrip/';
path(path, fieldtripfolder)
ft_defaults;
path('/home/kolai/Documents/Shishkina/ProgramFiles/Fieldtrip/external/mne/', path);


% This is an external matlab package used to save figures to PPTX
ft_defaults;
path('home/kolai/Documents/Shishkina/External_matlab/exportToPPTX_master/', path);

% EEEG lab functions
path (path, '/home/kolai/Documents/Shishkina/External_matlab/eeglab2019_0/functions/sigprocfunc');
path(path,'/home/kolai/Documents/Shishkina/External_matlab/eeglab2019_0/functions/popfunc');
path(path,'/home/kolai/Documents/Shishkina/External_matlab/eeglab2019_0/functions/guifunc');
path(path,'/home/kolai/Documents/Shishkina/External_matlab/eeglab2019_0/functions/adminfunc');
file = '/home/kolai/Documents/Shishkina/NeuralDataAnalysis/Autism/Code/scripts/for_meg_locs1.set'; % loc file
EEG = pop_loadset(file);
realdatapath = '/home/kolai/Documents/Shishkina/NeuralDataAnalysis/Autism/0101/MEGdata/';
DATAPATH  = '/home/kolai/Documents/Shishkina/NeuralDataAnalysis/Autism/0101/MEGdata/FT_beamf/';


%%
                              
%SUBJ = [ 'G008';  'G008'; 'G009'; 'G009'; 'G001'; 'G001';  'G002'; 'G003'; 'G003'; 'G005';'G006';'G007'; 'G007';  'G013'; 'G014'; 'G015'; 'G015'; 'G016'; 'G016'; 'G017'; 'G017'; 'G020'; 'G020'; 'G021'; 'G021'; 'G022'; 'G022'; 'G023'; 'G023'; 'G025'; 'G025'; 'G026'; 'G026'; 'G027'; 'G027'; 'G028'; 'G028'; 'G030'; 'G030'];
%DATE = ['f'; 'l';'f'; 'l'; 'f';    'l';      'l';       'f';    'l';    'l';   'f';   'f';   'l'; 'l'; 'l'; 'f'; 'l'; 'f'; 'l'; 'f'; 'l';  'f'; 'l'; 'f'; 'l'; 'f'; 'l'; 'f'; 'l';'f'; 'l'; 'f'; 'l'; 'f'; 'l'; 'f'; 'l'; 'f'; 'l'];
SUBJ = ['0101'];
DATE = ['oct-6'];
%%
s=1; 
size (SUBJ,1)
close all
subj = SUBJ (s,:); date = DATE (s,:);
DATAfolder = strcat(realdatapath);
savemegto = strcat(DATAPATH, subj, '_', date, '/', megfolder);
mkdir(savemegto)

%start PPTX report
exportToPPTX('new');
PPTXname  = strcat(savemegto, '/', subj, '_sensors_report');

    
    %%
SUBJ2 = ['0102'];
subj2 = SUBJ2(s,:)

%
ep_fiff_file = strcat(DATAfolder, subj2, '-noerror-lagcorrected-epo.fif');
hdr = ft_read_header(ep_fiff_file);
ev_name=[realdatapath, 'epochs1/', subj2, date, '_events.mat']
load(ev_name)
% load epochs
cfg = [];
cfg.dataset = ep_fiff_file;
cfg.channel={'MEG*'};
epochs = ft_preprocessing(cfg);
%  select epochs according to events

ev1 = find(events(:,3)==2);
ev2 = find(events(:,3)==4);
ev3 = find(events(:,3)==8);
EV = {ev1, ev2, ev3};
    
     %% Load raw data
%     fiff_file = strcat(realdatapath, subj, '/ICA_notch_', date, '/', subj, '_static_raw.fif');
%     hdrraw = ft_read_header(fiff_file);
%     first= round(cast(hdrraw.orig.raw.first_samp, 'double'));
%     events(:,1) =  events(:,1)-first;
%     
%     % Define trials
%     % trl:   start, end and offset (interval before the event) 
%     pre = -1.0* hdr.Fs ;
%     post = 1.2* hdr.Fs ;
%     trl=[];
%     for i=1:size (events,1)
%         trl(i, 1)=(events(i,1)+pre) ; 
%         trl(i, 2)=(events(i,1)+post) ; 
%         trl(i, 3)= -1.0*hdr.Fs ; % offset
%         trl(i, 4) = events(i,3); % stimulus_value;
%     end
% 
%     % extract data and epochs from the raw
%     cfg = [];
%     cfg.trl=trl;
%     cfg.channel     = 'meg';
%     % cfg.demean      = 'yes';
%     cfg.dftfilter   = 'yes';
%     cfg.dftfreq     = [50 100];
%     cfg.dataset = fiff_file;
%     cfg    = ft_definetrial(cfg);
%     epochs = ft_preprocessing(cfg);
%     
%% Browse epochs
% %     cfg=[]
% %     cfg.viewmode = 'vertical';
% %     cfg.plotevents = 'yes';
% %     cfg.blocksize =5;
% %     ft_databrowser(cfg, epochs);
    
%% Plot average ERF
cfg = [];
cfg.channel=epochs.label;
avg = ft_timelockanalysis(cfg,epochs);
%%
hh=figure;
ttt = find(avg.time>-0.3 & avg.time<0.5);
plot (avg.time(ttt), avg.avg(:, ttt));
title ('Sensor averages')
exportToPPTX('addslide'); % slide with distributions
exportToPPTX('addpicture', hh, 'Position', [0.5,0.5,9,6]);
close (hh)

%% Calculating power in sensors, find max weighted freq, freq range and P of gamma enhancement in this freq range.
hh=figure;
LineColor{1} = 'k'; LineColor{2} = 'b'; LineColor{3} = 'g'; 
LineColorMax{1} = '--k'; LineColorMax{2} = '--b'; LineColorMax{3} = '--g'; 
ax=gca;

Ch = {'MEG1932',  'MEG1922', 'MEG2042',  'MEG2032',  'MEG2112', 'MEG2122',  'MEG2342', 'MEG2332',  'MEG1732', 'MEG1942', 'MEG1912', 'MEG2012', 'MEG2022', 'MEG2312', 'MEG2322', 'MEG2512',...
      'MEG1933',  'MEG1923', 'MEG2043',  'MEG2033',  'MEG2113', 'MEG2123',  'MEG2343', 'MEG2333',  'MEG1733', 'MEG1943', 'MEG1913', 'MEG2013', 'MEG2023', 'MEG2313', 'MEG2323', 'MEG2513'};

exportToPPTX('addslide'); % slide with distributions
exportToPPTX('addtext','Power ratio [post/pre] distribution for GRA1 and GRA2, conditions 1,2,3');
exportToPPTX('addtext','GRA1 sensor set, ... MEG0112,...', 'Position', [0.5, 0.5, 10, 1]);
exportToPPTX('addtext','GRA2 sensor set, ... MEG0113,...', 'Position', [0.5, 3.6, 10, 1]);

X=[0.4, 2.6, 4.8, 7.0; 0.4, 2.6, 4.8, 7.0]; Y=[1, 1, 1, 1; 4, 4, 4, 4]; % positions for plot

    %%
for con=1:3 % for conditions
    cfg = [];
    cfg.trials = EV{con}; % EV{con}';
    [selepo] = ft_selectdata(cfg, epochs);

    cfg = [];
    cfg.channel=epochs.label;
    avg = ft_timelockanalysis(cfg,selepo);

   % subtract evoked! Important for Moscow data with 60Hz projector.
    for ttt=1:size (selepo.trial,2)
        selepo.trial{ttt}=selepo.trial{ttt}-avg.avg;
    end     

    cfg.latency = [-0.9 0.0];
    latency_pre = cfg.latency;
    [dataPre] = ft_selectdata(cfg, selepo);
    cfg.latency = [0.3 1.2];
    latency_post = cfg.latency;
    [dataPost] = ft_selectdata(cfg, selepo);

    cfg = [];
    cfg.method    = 'mtmfft';
    cfg.output    ='pow'; % 'fourier'
    cfg.taper        = 'dpss';
    cfg.keeptrials = 'yes';
    cfg.tapsmofrq = tapsmofrq;
    cfg.foilim    = gammarange; %freq band of interest
    freqPre = ft_freqanalysis(cfg, dataPre);  % trials x Ch x freq
    freqPost = ft_freqanalysis(cfg, dataPost);
    % N of max channel
    for j=1:length(Ch)
        [ch(j),x] = find(strcmp(freqPost.label,Ch{j}));
    end
    %% Find peak freq based on max probability value
    % b1 is the max ch
    [a,b1] = max(  squeeze(mean(  mean(freqPost.powspctrm(:,ch,:),3),1))./squeeze(mean(mean(freqPre.powspctrm(:,ch,:),3),1))  ); % trial x ch x freq
    ratio{con} =  (  squeeze(mean(freqPost.powspctrm(:,ch(b1),:),1))./squeeze(mean(freqPre.powspctrm(:,ch(b1),:),1) ) -1   )'; % post/pre ratio at the 'max channel'
    ChName{con}  = freqPost.label(ch(b1));

    % find N=4 max channels
    change =squeeze(mean(freqPost.powspctrm(:,ch,:),1))./squeeze(mean(freqPre.powspctrm(:,ch,:),1)); % ch x freq;
    change = mean(change,2); % over frequencies in gamma range
    change =[1:length(change)'; change']';
    change = sortrows(change,2);
    selection4=change([end, end-1, end-2, end-3],1); %N channels with max (max-1, max-2, etc) power CHANGE in the GAMMA BAND
    t= ch (selection4);

    Ratio{con} =  (  squeeze(mean(mean(freqPost.powspctrm(:,t,:),1),2))./squeeze(mean(mean(freqPre.powspctrm(:,t,:),1),2) ) -1   )'; % post/pre ratio at the 'max channel selection'
    plot (ax,   freqPost.freq, Ratio{con}, LineColor{con}); hold on, plot(ax, freqPost.freq, ratio{con},  LineColorMax{con});

    bins{con} = [];
    for f=1:size(freqPost.powspctrm,3) % for all freq
         post = squeeze(freqPost.powspctrm(:,ch(b1),f)); % ch(b1) is the max ch
         pre  = squeeze(freqPre.powspctrm(:,ch(b1),f));    
         [p_f{con}(f),h,stats] = ranksum(pre,post);
    end
    [val, ind ] = min(p_f{con})
    PeakF_Pbased{con}(1) = freqPre.freq(ind);
    PeakF_Pbased{con}(2) = val;

    %% Find peak freq based on absolute max value

    %find max weighted freq
    [MAXratio, MAXind] = (max(Ratio{con})); maxF = freqPost.freq(MAXind);
    bins{con} = find( (Ratio{con}) > max(Ratio{con} )/3*2); % 'bins' is the subject/condition-specific band
    binsfreqs{con} = freqPost.freq(bins{con});

    % if frequencies are too far from the maximum (more than 20 Hz) we exclude them!
    excl = find(abs(binsfreqs{con}-maxF)>20);
    bins{con}(excl) =[];
    binsfreqs{con} = freqPost.freq(bins{con});

    wF{con} = sum(binsfreqs{con}.*Ratio{con}(bins{con}))/sum(Ratio{con}(bins{con}));  % if 'Ratio', then at the 'max channel selection, N=4, see line 171'
    wPmax{con} = mean(ratio{con}(bins{con})+1); % max ch
    wPsel{con} = mean(Ratio{con}(bins{con})+1); % selection of 4 channels

    % in the max channels find P of gamma increase
    [aa,bb] = min(abs(freqPost.freq-wF{con})); % bb is the closest max freq    

    % for subject/condition-specific band
    post = squeeze(mean(mean(freqPost.powspctrm(:,t,bins{con}),3),2));  %bins{con} t is the max ch selection (of 4 ch), at the subject/condition-specific band
    pre  = squeeze(mean(mean(freqPre.powspctrm(:,t,bins{con}),3) ,2));
    [p_band{con},h,stats] = ranksum(pre,post);
    postave{con} = mean(post)*1e24;
    preave{con}  = mean(pre)*1e24;

    % for maximum
    post = squeeze(mean(freqPost.powspctrm(:,t,bb),2));  %bb is the nearest to wF frequency bin
    pre = squeeze(mean(freqPre.powspctrm(:,t,bb),2));       
    [p_max{con},h,stats] = ranksum(pre,post);

%       postave{con} = squeeze(mean(mean(freqPost.powspctrm(:,(ch(b)),bins{con}),3),1))*1e24;
%       preave{con}  = squeeze(mean(mean(freqPre.powspctrm(:,(ch(b)),bins{con}),3),1))*1e24;
%%
    data = squeeze(mean(mean(freqPost.powspctrm(:,:,ind),3),1))./squeeze(mean(mean(freqPre.powspctrm(:,:,ind),3),1))-1; % ind is the p-based peak freq
    [gra,y] = find(strcmp(freqPost.grad.chantype, 'megplanar'));

    data = data(gra);
    Data=data(1:2:length(data));
    h=figure;
    topoplot (Data, EEG.chanlocs, 'electrodes', 'on', 'maplimits', [0,max(ratio{1})]);    colorbar
    title (strcat(subj, ',..GRA1..V', num2str(con), ',..maxCh=', ChName{con}),  'FontSize', 20)
    exportToPPTX('addpicture', h, 'Position', [X(1,con),Y(1,con),2.3,1.8]);
    Data=data(2:2:length(data));
    close (h)
    h=figure;
    topoplot (Data, EEG.chanlocs, 'electrodes', 'on','maplimits', [0,max(ratio{1})]);    colorbar
    title (strcat(subj, ',..GRA2..V', num2str(con), ',..maxCh=', ChName{con}), 'FontSize', 20)
    exportToPPTX('addpicture', h, 'Position', [X(2,con),Y(2,con),2.3,1.8]);
    close (h)
end  % end for conditions in sensor space
%%
title(ax, strcat('[post-pre]/pre ratio in MAX-prob and 3-AVE channels:',  ChName{1}, '../',ChName{2}, '../', ChName{3}),  'FontSize', 10 )
legend(ax, {'SlowAVE', 'SlowMAX', 'Medium','MediumMAX',  'FastAVE', 'FastMAX'}, 'FontSize', 10);
exportToPPTX('addslide'); % slide with distributions
exportToPPTX('addpicture', hh, 'Position', [1,1,6,4]);
exportToPPTX('addtext', strcat('F1=', num2str(wF{1}), ', F2=', num2str(wF{2}), ', F3=', num2str(wF{3})), 'FontSize', 11, 'Position', [0.5, 5.0, 10, 1]);
exportToPPTX('addtext', strcat('P1max_dB=', num2str(wPmax{1}),'P2max_dB=', num2str(wPmax{2}), ', P3max_dB=', num2str(wPmax{3})), 'FontSize', 11, 'Position', [0.5, 5.4, 10, 1]);
exportToPPTX('addtext', strcat('P1selection4_dB=', num2str(wPsel{1}), ', P2selection4_dB=', num2str(wPsel{2}), ', P3selection4_dB=', num2str(wPsel{3})), 'FontSize', 11, 'Position', [0.5, 5.6, 10, 1]);

exportToPPTX('addtext', strcat('Condition-specific band for V1:...', num2str(binsfreqs{1}(1)), '..to..', num2str(binsfreqs{1}(length(binsfreqs{1})))   ), 'FontSize', 11,  'Position', [0.5, 6.0, 10, 1]);
exportToPPTX('addtext', strcat('Condition-specific band for V2:...', num2str(binsfreqs{2}(1)), '..to..', num2str(binsfreqs{2}(length(binsfreqs{2})))   ), 'FontSize', 11,  'Position', [0.5, 6.2, 10, 1]);
exportToPPTX('addtext', strcat('Condition-specific band for V3:...', num2str(binsfreqs{3}(1)), '..to..', num2str(binsfreqs{3}(length(binsfreqs{3})))   ), 'FontSize', 11,  'Position', [0.5, 6.4, 10, 1]);

exportToPPTX('addtext', strcat('M-W rank test at freq closest to wF:  pV1_max=', num2str(p_max{1}), ', pV2_max=', num2str(p_max{2}), ', pV3_max=', num2str(p_max{3})), 'FontSize', 11,  'Position', [0.5, 7.1, 10,1]);
exportToPPTX('addtext', strcat('M-W rank test in subj/freeq-specific band:  pV1_band=', num2str(p_band{1}), ', pV2_band=', num2str(p_band{2}), ', pV3_band=', num2str(p_band{3})), 'FontSize', 11, 'Position', [0.5, 7.3, 10,1]);
close (gcf)


%% Save sensor info
readme = [];
readme.wF = 'Individual weighterd frequency, at 2/3 max in max ch';
readme.wPmax = 'Max weighted pow in dB in max channel ';
readme.wPsel = 'Max weighted pow in dB in the selection-of-4-max-channels';
readme.p_band = 'probability at the selection-of-4-max-channels amplitude increase ';
readme.p_max = 'probability at the selection-of-4-max-channels amplitude increase, is the nearest to wF frequency bin ';
readme.ChName = 'Max sensor name';
readme.preave = 'prestim power in max sensor selection, averaged over condition-specific band';
readme.postave = 'poststim power in max sensor selection, averaged over condition-specific band';
readme.binsfreqs = 'condition specific bands';
readme.PeakF_Pbased= 'peak F at max probability and this probability';

Savename = strcat(savemegto,'/', subj, '_sensors_res.mat');
% Max weighted freq
save (Savename, 'wF');
% Max weighted pow in dB in max channel
save (Savename, 'wPmax', '-append');
% Max weighted pow in dB in the channel selection (n=4)
save (Savename, 'wPsel', '-append');
% P value for condition-specific band
save (Savename, 'p_band', '-append');
% P value for the nearest max freq
save (Savename, 'p_max', '-append');
% Max Sensor
save (Savename, 'ChName', '-append');
% prestim power in max sensor selection, averaged over condition-specific band
save (Savename, 'preave', '-append'); 
% poststim power in max sensor selection, averaged over condition-specific band
save (Savename, 'postave', '-append'); 
% condition specific bands
save (Savename, 'binsfreqs', '-append');   
save (Savename, 'PeakF_Pbased', '-append');   
save (Savename, 'latency_post', 'latency_pre', '-append');   

%%
exportToPPTX('saveandclose', PPTXname);




