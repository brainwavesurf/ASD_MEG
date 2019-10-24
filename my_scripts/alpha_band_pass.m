%% to apply an alpha bandpass filter to the whole continuous recording 
%%
clear all;
close all;

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(fieldtripfolder, path);
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);
realdatapath = '/home/a_shishkina/data/KI/SUBJECTS/';
savepath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';

%%
%add list of subjects:
SUBJ = [ '0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0135'; '0136'; '0137'; '0138'; '0139';...  
         '0140'; '0141'; '0158'; '0159'; '0160'; '0161'; '0162'; '0163'; '0164'; '0178'; '0179'; '0253';...
         '0254'; '0255'; '0256'; '0257'; '0259'; '0273'; '0274'; '0275'; '0276'; '0277'; '0278'; '0278';...
         '0310'; '0346'; '0347'; '0348'; '0349'; '0350'; '0351'; '0357'; '0358'; '0378'; '0379'; '0380'; '0381';...
         '0382'; '0383'; '0384']; 
%% loop for all subjects
for s=1: size (SUBJ,1)
    
    subj = SUBJ (s,:); 
    
    rawfile = strcat(realdatapath, subj, '/ICA_nonotch_crop/', subj, '_rings_ICA_raw.fif');
    hdrraw = ft_read_header(rawfile);
    
    cfg =[];
    cfg.dataset = rawfile;
    cfg.bpfreq = [8 13];
    [alpha_bp] = ft_preprocessing(cfg);
    
    %save alpha band-passed recordings
    filename = strcat(realdatapath, subj, '/ICA_nonotch_crop/', subj, '_alpha_band_pass.mat');
    save(filename, 'alpha_bp', '-v7.3');
end
