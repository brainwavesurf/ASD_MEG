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

%SUBJ = ['0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0136'; '0137'; '0138'; '0139'; '0140'; '0141'; '0158'; '0159'; '0160'; '0161'; '0162'; '0163'; '0164'; '0178'; '0179'; '0253'; '0254'; '0255'; '0256'; '0257'; '0259'; '0273'; '0274'; '0275'; '0276'; '0277'; '0346'; '0347'; '0348'; '0350'; '0351'; '0357'; '0358'; '0378'; '0380'; '0381'; '0382'; '0383'; '0384'; '0385']; 
SUBJ = ['0076'];

%%
for s=1: size (SUBJ,1)
    close all
    subj = SUBJ (s,:);
    mkdir (strcat(savedatapath , subj));
    savefolder = (strcat(savedatapath , subj));

    %% LOAD EVENTS: clean events order number
    load([visdatapath, subj, '/ICA_nonotch_crop/', subj, '_clean_events.mat']) 
    
    %% for what?
    if events(length(events),1)/SF/60<19 || events(length(events),1)/SF/60>35
        'WRONG SF!!!' 
        return
    end
    %%  select visual epochs according to events
    ev1 = find(events(:,3)==2);
    ev2 = find(events(:,3)==4);
    ev3 = find(events(:,3)==8);
    EV = {ev1, ev2, ev3}; 
    
    info = [];
    for eve=1:3
        info.evnum{eve} = clean_events_N(EV{eve}); %clean_events_N - order number of the event in original dataset (read from .mat file above)
        info.rt{eve} = RT_good(EV{eve});
        prev_stim_type=[]; 
        prev_stim_length=[];
        for i = 1:size(EV{eve},1)
            evnum = clean_events_N(EV{eve}(i));
            if (evnum~=1) && (evnum~=91) && (evnum~=181) && (evnum~=226) %if not the 1st event of the run
                prev_stim_length(i) = Stim_length(evnum-1); 
                prev_stim_type(i) = eventtype(evnum-1);
            else
                prev_stim_length(i) = 0;
                prev_stim_type(i) = 0;
            end
        end
        info.prev_stim_type{eve} = prev_stim_type;
        info.prev_stim_length{eve} = prev_stim_length;
    end 

    %% The same without dividing to stim types:
    allinfo=[];
    
    for i=1: size(events,1)
       allinfo.event(i) = events(i,3); 
       allinfo.evnum(i) = clean_events_N(i);
       allinfo.rt(i)=RT_good(i);
       evnum=clean_events_N(i);
       if (evnum~=1) && (evnum~=91)  && (evnum~=181) && (evnum~=226)
           allinfo.prev_stim_type(i)=eventtype(evnum-1);
           allinfo.prev_stim_length(i)=Stim_length(evnum-1);
       end
    end

    %% save stim info
    savefile = [savefolder, '/',  subj, '_info'];
    save (savefile , 'info', 'allinfo');

end % for subjects
%%
   




 