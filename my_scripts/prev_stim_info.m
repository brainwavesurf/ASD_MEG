% step0_stim_info.m
% Keep info about preceding stimulus, write structures to matspaces, individually for each subject:
% allinfo - info about all good trials
% info - info about good trial organized according to the stim type (1-slow/2-medium/3-fast)
%%
clear
close all
screensize = get( groot, 'Screensize' );

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(path, fieldtripfolder)
ft_defaults;

visdatapath = '/home/a_shishkina/data/KI/SUBJECTS/';
savedatapath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';

SUBJ = ['0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0136'; '0137'; '0138'; '0139'; '0140'; '0141'; '0158'; '0159'; '0160'; '0161'; '0162'; '0163'; '0164'; '0178'; '0179'; '0253'; '0254'; '0255'; '0256'; '0257'; '0259'; '0273'; '0274'; '0275'; '0276'; '0277'; '0346'; '0350'; '0351'; '0357'; '0358'; '0378'; '0380'; '0381'; '0382'; '0383'; '0384']; 
%SUBJ = ['0076'];
%'0347';'0348';
%%
for s=1: size (SUBJ,1)
    close all
    subj = SUBJ (s,:);
    mkdir (strcat(savedatapath , subj));
    savefolder = (strcat(savedatapath , subj));
    %% load all events
    all = load([visdatapath, subj, '/ICA_nonotch_crop/', subj,'_events.mat']);
    
    %% load clean events
    load([visdatapath, subj, '/ICA_nonotch_crop/', subj, '_clean_events.mat']); 
    
    %% time to previous event in sec
    for i = 2:size(all.events,1) % 270 events
        all_to_prev(i) = double(all.events(i,1) - all.events(i-1,1))/SF;
    end
    
    %% allinfo contains info about good events
    allinfo = [];
    for i = 1: size(events,1) % 227 good events
        allinfo.event(i) = events(i,3); % event in this trial
        allinfo.epo_num(i) = clean_events_N(i);
        allinfo.rt(i) = RT_good(i);
        epo_num = clean_events_N(i);
        if (epo_num~=1) && (epo_num~=91)  && (epo_num~=181) && (epo_num~=226) && (all_to_prev(epo_num)<5.3) % in sec, should be shorter then 3+1+0.1+1.2
           allinfo.prev_stim_type(i) = eventtype(epo_num-1); % event in previous trial
           allinfo.prev_stim_length(i) = Stim_length(epo_num-1);
           allinfo.time_to_previous_ev(i) = all_to_prev(epo_num);
        end
    end
    %% save stim info
    savefile = [savefolder, '/',  subj, '_info'];
    save (savefile , 'allinfo');

end % for subjects
%%
   




 