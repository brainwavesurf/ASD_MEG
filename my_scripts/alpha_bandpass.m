%% Apply band-pass filter to the whole continious recording
% Cut epochs [-0.8 0]

clear all;
close all;

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(fieldtripfolder, path);
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);

realdatapath = '/home/a_shishkina/data/KI/SUBJECTS/';
savepath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';

%%
% add list of all subjects:
SUBJ_All = ['0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0135'; '0136'; '0137'; '0138'; '0139';...  
            '0140'; '0141'; '0158'; '0159'; '0160'; '0161'; '0162'; '0163'; '0164'; '0178'; '0179'; '0253';...
            '0254'; '0255'; '0256'; '0257'; '0259'; '0273'; '0274'; '0275'; '0276'; '0277'; '0278'; '0278';...
            '0310'; '0346'; '0347'; '0348'; '0349'; '0350'; '0351'; '0357'; '0358'; '0378'; '0379'; '0380'; '0381';...
            '0382'; '0383'; '0384'];

% choose only subjects with MRI 
     
SUBJ_NT = [ '0101'; '0102'; '0103'; '0104'; '0105'; '0135'; '0136';...  
            '0137'; '0138'; '0140'; '0158'; '0162'; '0163'; '0178';...
            '0179'; '0255'; '0257'; '0348'; '0378'; '0379'; '0384']; 
        
SUBJ_ASD = ['0106'; '0107'; '0139'; '0141'; '0159'; '0160'; '0161';...  
            '0164'; '0253'; '0254'; '0256'; '0273'; '0274'; '0275';...
            '0276'; '0346'; '0347'; '0349'; '0351'; '0358';...
            '0380'; '0381'; '0382'; '0383'];  
%without '0357';

SUBJ = [SUBJ_NT; SUBJ_ASD];

for s=1: size (SUBJ,1)
    close all
    subj = SUBJ (s,:); 
    savemegto = strcat(savepath, subj);
    epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
 
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

    % extract epochs from alpha band-passed raw data
    cfg = [];
    cfg.trl         = trl;
    cfg.dataset     = fiff_file;
    cfg.channel     = 'MEG';
    cfg.bpfilter    = 'yes';
    cfg.bpfreq      = [10 14];
    %cfg.bpfilttype = 'firws';
    
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
    slow_alpha_bp = ft_selectdata(cfg, epochs); %extraxt trials after slow stimuli
    
    cfg.trials = medium_ind;
    cfg.latency = [-0.8 0];
    medium_alpha_bp = ft_selectdata(cfg, epochs); %extraxt trials after slow stimuli
    
    cfg.trials = fast_ind;
    cfg.latency = [-0.8 0];
    fast_alpha_bp = ft_selectdata(cfg, epochs); %extraxt trials after fast stimuli
    
    filename = strcat(epofolder, subj, '_preproc_alpha_bp_epochs.mat');
    save(filename, 'epochs', 'slow_alpha_bp', 'medium_alpha_bp', 'fast_alpha_bp');
    
end
