%% to apply an hp = 35 filter to the whole continuous recording and cut for epochs according slow/medium/fast condition 
%and interstimuli and stimulation interval
%%
clear all;
close all;

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(fieldtripfolder, path);
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
%% loop for all subjects
for s=1: size (SUBJ,1)
    
    subj = SUBJ (s,:); 
    savemegto = strcat(savepath, subj);
    epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
    
    %% load group preceding events in order to select the epochs later on
    load ([ savemegto, '/', subj, '_info.mat']);
    
    slow_ind = find(allinfo.prev_stim_type==2);
    medium_ind = find(allinfo.prev_stim_type==4);
    fast_ind = find(allinfo.prev_stim_type==8);

  %% Load unfiltered epochs and divide them according to preceding conditions
    
    % load raw fif data after ICA 
    fiff_file = strcat(realdatapath, subj, '/ICA_nonotch_crop/', subj, '_rings_ICA_transstandard-raw.fif');
    hdrraw = ft_read_header(fiff_file);
    
    % load info about events
    load(strcat(realdatapath, subj, '/ICA_nonotch_crop/', subj, '_clean_events.mat'));
     
    first = round(cast(hdrraw.orig.raw.first_samp, 'double'));
    events(:,1) =  events(:,1) - first;
    
    %Define trials
    % trl:   start, end and offset (interval before the event)
    pre = -1.0* hdrraw.Fs ;
    post = 1.2* hdrraw.Fs ;
    trl=[];
    for i=1:size (events,1)
        trl(i, 1)=(events(i,1)+pre) ;
        trl(i, 2)=(events(i,1)+post) ;
        trl(i, 3)= -1.0*hdrraw.Fs ; % offset
        trl(i, 4) = events(i,3); % stimulus_value;
    end
    
    filename = strcat(epofolder, subj, '_trials.mat');
    save(filename, 'trl');
    
    % extract epochs from alpha band-passed raw data
    cfg = [];
    cfg.trl         = trl;
    cfg.dataset     = fiff_file;
    cfg.channel     = 'MEG';
    cfg.hpfreq      = 35;
    cfg.demean      = 'yes';
    
    [cfg] = ft_definetrial(cfg);
    epochs = ft_preprocessing(cfg);
    
    %% load info about preceding events in order to select the epochs later on
    load ([ savemegto, '/', subj, '_info.mat'])
    
    slow_ind = find(allinfo.prev_stim_type==2); %slow
    medium_ind = find(allinfo.prev_stim_type==4); %medium
    fast_ind = find(allinfo.prev_stim_type==8); %fast
    
    cfg = [];
    cfg.trials = slow_ind;
    cfg.latency = [-0.8 0];
    slow_alpha_pre = ft_selectdata(cfg, epochs); %extraxt trials after slow stimuli
    cfg.latency = [0.4 1.2];
    slow_alpha_post = ft_selectdata(cfg, epochs);
    
    cfg.trials = medium_ind;
    cfg.latency = [-0.8 0];
    medium_alpha_pre = ft_selectdata(cfg, epochs); %extraxt trials after slow stimuli
    cfg.latency = [0.4 1.2];
    medium_alpha_post = ft_selectdata(cfg, epochs); %extraxt trials after slow stimuli
    
    cfg.trials = fast_ind;
    cfg.latency = [-0.8 0];
    fast_alpha_pre = ft_selectdata(cfg, epochs); %extraxt trials after fast stimuli
    cfg.latency = [0.4 1.2];
    fast_alpha_post = ft_selectdata(cfg, epochs); %extraxt trials after fast stimuli
    
    filename = strcat(epofolder, subj, '_alpha_epochs.mat');
    save(filename, 'epochs', 'epochs_pre', 'epochs_post', 'slow_alpha_pre', 'slow_alpha_post', 'medium_alpha_pre', 'medium_alpha_post', 'fast_alpha_pre', 'fast_alpha_post');  
end
