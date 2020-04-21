%% to apply a band-pass = 2-40 Hz filter to the whole continuous recording 
%to cut epochs according slow/medium/fast condition 
%in the interstimulus [-0.8 0] and the stimulation [0.4 1.2] intervals
%%
close all;

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(fieldtripfolder, path);
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);
pathfrom = '/net/server/data/Archive/aut_gamma/orekhova/KI/';
realdatapath = '/net/server/data/Archive/aut_gamma/orekhova/KI/SUBJECTS/';
savepath = '/net/server/data/Archive/aut_gamma/orekhova/KI/Scripts_bkp/Shishkina/KI/Results_Alpha_and_Gamma/';

%%
%load subj list
     
SUBJ_NT = [ '0101'; '0102'; '0103'; '0104'; '0105'; '0135'; '0136';...  
            '0137'; '0138'; '0140'; '0158'; '0162'; '0163'; '0178';...
            '0179'; '0255'; '0257'; '0348'; '0378'; '0379'; '0384']; 
        
SUBJ_ASD = ['0106'; '0107'; '0139'; '0141'; '0159'; '0160'; '0161';...  
            '0164'; '0253'; '0254'; '0256'; '0273'; '0274'; '0275';...
            '0276'; '0346'; '0347'; '0349'; '0351'; '0358';...
            '0380'; '0381'; '0382'; '0383'];  

SUBJ = [SUBJ_ASD; SUBJ_NT];

for s=1:size (SUBJ,1)

    subj = SUBJ (s,:); 
    
    % load raw fif data after ICA 
    fiff_file = strcat(realdatapath, subj, '/ICA_nonotch_crop/', subj, '_rings_ICA_raw.fif');
    hdrraw = ft_read_header(fiff_file);
    
    % load info about events
    load(strcat(realdatapath, subj, '/ICA_nonotch_crop/', subj, '_clean_events.mat'));
     
    first = round(cast(hdrraw.orig.raw.first_samp, 'double'));
    events(:,1) =  events(:,1) - first;
    
    %Define trials
    % trl:   start, end and offset (interval before the event)
    pre = -1.0* hdrraw.Fs ;
    post = 1.2* hdrraw.Fs ;
    trl = zeros(size(events,1),4);
    for i=1:size(events,1)
        trl(i, 1) = (events(i,1)+pre) ;
        trl(i, 2) = (events(i,1)+post) ;
        trl(i, 3) = -1.0*hdrraw.Fs ; % offset
        trl(i, 4) = events(i,3); % stimulus_value;
    end

    % extract epochs from alpha band-passed raw data
    cfg = [];
    cfg.trl         = trl;
    cfg.dataset     = fiff_file;
    cfg.channel     = 'MEG';
    cfg.bpfilter    = 'yes';
    cfg.bpfreq      = [2 40];
    
    [cfg] = ft_definetrial(cfg);
    epochs = ft_preprocessing(cfg);
    
    %% load info about preceding events in order to select the epochs later on
    load(strcat(pathfrom, '/Results_Alpha_and_Gamma/', subj, '/',  subj, '_info.mat'))
    
    slow_ind = find(ALLINFO.stim_type_previous_tr==2); %slow
    medium_ind = find(ALLINFO.stim_type_previous_tr==4); %medium
    fast_ind = find(ALLINFO.stim_type_previous_tr==8); %fast
    
    %select epochs in different conditions
    cfg = [];
    cfg.trials = slow_ind;
    slow_alpha_epochs = ft_selectdata(cfg, epochs); 
    
    cfg.trials = medium_ind;
    medium_alpha_epochs = ft_selectdata(cfg, epochs);  
    
    cfg.trials = fast_ind;
    fast_alpha_epochs = ft_selectdata(cfg, epochs);
    
    %select interstimulus and stimulation intervals
    cfg = [];
    cfg.trials = slow_ind;
    cfg.latency = [-0.8 0];
    slow_alpha_isi = ft_selectdata(cfg, epochs); %extraxt trials after slow stimuli
    
    cfg.trials = medium_ind;
    cfg.latency = [-0.8 0];
    medium_alpha_isi = ft_selectdata(cfg, epochs); %extraxt trials after slow stimuli
    
    cfg.trials = fast_ind;
    cfg.latency = [-0.8 0];
    fast_alpha_isi = ft_selectdata(cfg, epochs); %extraxt trials after fast stimuli
    
    cfg = [];
    cfg.trials = slow_ind;
    cfg.latency = [0.4 1.2];
    slow_alpha_post = ft_selectdata(cfg, epochs); %extraxt trials after slow stimuli
    
    cfg.trials = medium_ind;
    cfg.latency = [0.4 1.2];
    medium_alpha_post = ft_selectdata(cfg, epochs); %extraxt trials after slow stimuli
    
    cfg.trials = fast_ind;
    cfg.latency = [0.4 1.2];
    fast_alpha_post = ft_selectdata(cfg, epochs);
     
    filename = strcat(savepath, subj, '/', subj, '_preproc_alpha_2_40_epochs.mat');
    save(filename, 'epochs', 'slow_alpha_epochs','medium_alpha_epochs','fast_alpha_epochs', ...
        'slow_alpha_isi', 'medium_alpha_isi', 'fast_alpha_isi', 'slow_alpha_post', 'medium_alpha_post', 'fast_alpha_post');
end
