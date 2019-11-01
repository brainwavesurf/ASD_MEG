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
SUBJ = ['0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0136'; '0137'; '0138'; '0139'; '0140'; '0141'; '0158'; '0159'; '0160'; '0161'; '0162'; '0163'; '0164'; '0178'; '0179'; '0253'; '0254'; '0255'; '0256'; '0257'; '0259'; '0273'; '0274'; '0275'; '0276'; '0277'; '0346'; '0347'; '0348'; '0350'; '0351'; '0357'; '0358'; '0378'; '0380'; '0381'; '0382'; '0383'; '0384']; 


for s=1: size (SUBJ,1)
    close all
    subj = SUBJ (s,:); 
    savemegto = strcat(savepath, subj);
    epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
    
    %% load group preceding events in order to select the epochs later on
    load ([ savemegto, '/', subj, '_info.mat'])
    
    slow_ind = find(allinfo.prev_stim_type==2);
    medium_ind = find(allinfo.prev_stim_type==4);
    fast_ind = find(allinfo.prev_stim_type==8);
 
    %% Load unfiltered epochs and divide them according to preceding conditions
    ep_fiff_file = strcat(epofolder, subj, '-noerror-lagcorrected-epo.fif')
    hdr = ft_read_header(ep_fiff_file);
    
    cfg           = [];  
    cfg.dataset   = ep_fiff_file;
    cfg.channel   = {'MEG'};
    epochs        = ft_preprocessing(cfg);
    
    cfg           = [];
    cfg.trials    = slow_ind;
    slow_epochs   = ft_selectdata(cfg, epochs); %extraxt trials after slow stimuli
    
    cfg.trials    = fast_ind;
    fast_epochs   = ft_selectdata(cfg, epochs); %extraxt trials after fast stimuli
    
    filename = strcat(epofolder, subj, '_preproc_epochs.mat');
    save (filename, 'epochs', 'slow_epochs', 'fast_epochs');
    
end