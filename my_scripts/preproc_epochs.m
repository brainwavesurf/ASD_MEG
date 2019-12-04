%%
clear all;
close all;

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(path, fieldtripfolder)
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);

realdatapath = '/home/a_shishkina/data/KI/SUBJECTS/';
savepath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';

%add list of subjects:
SUBJ = ['0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0136'; '0137'; '0138'; '0139'; '0140'; '0141'; '0158'; '0159'; '0160'; '0161'; '0162'; '0163'; '0164'; '0178'; '0179'; '0253'; '0254'; '0255'; '0256'; '0257'; '0259'; '0273'; '0274'; '0275'; '0276'; '0277'; '0346'; '0347'; '0348'; '0350'; '0351'; '0358'; '0378'; '0380'; '0381'; '0382'; '0383'; '0384']; 

% without SUBJ = ['0357'];

for s=1: size (SUBJ,1)
    close all
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
    pre = -1.0* hdr.Fs ;
    post = 1.2* hdr.Fs ;
    trl=[];
    for i=1:size (events,1)
        trl(i, 1)=(events(i,1)+pre) ;
        trl(i, 2)=(events(i,1)+post) ;
        trl(i, 3)= -1.0*hdr.Fs ; % offset
        trl(i, 4) = events(i,3); % stimulus_value;
    end

    % extract data and epochs from the raw
    cfg = [];
    cfg.trl=trl;
    cfg.channel     = 'meg';
    % cfg.demean      = 'yes';
    cfg.dftfilter   = 'yes';
    cfg.dftfreq     = [50 100];
    cfg.dataset = fiff_file;
    cfg    = ft_definetrial(cfg);
    epochs = ft_preprocessing(cfg);
    
%check approximately 100 ms response
%     cfg = [];
%     cfg.channel=epochs.label;
%     avg = ft_timelockanalysis(cfg,epochs);
%     
%     hh=figure;
%     ttt = find(avg.time>-0.3 & avg.time<0.5);
%     plot (avg.time(ttt), avg.avg(:, ttt));
%     

    cfg           = [];
    cfg.trials    = slow_ind;
    slow_epochs   = ft_selectdata(cfg, epochs); %extraxt trials after slow stimuli
    
    cfg.trials    = fast_ind;
    fast_epochs   = ft_selectdata(cfg, epochs); %extraxt trials after fast stimuli
    
    filename = strcat(epofolder, subj, '_preproc_epochs.mat');
    save (filename, 'epochs', 'slow_epochs', 'fast_epochs');
    
end