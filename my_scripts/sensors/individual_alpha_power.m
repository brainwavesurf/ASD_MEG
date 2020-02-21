%calculation individual alpha freq with linear trend removal
clear all;
close all;

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(path, fieldtripfolder)
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);
path('/home/a_shishkina/projects/ASD_MEG/Git/ASD_MEG/my_scripts', path);
path('/home/a_shishkina/externals/', path);

savepath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';
realdatapath = '/home/a_shishkina/data/KI/SUBJECTS/';

%load subj list
SUBJ_NT = [ '0101'; '0102'; '0103'; '0104'; '0105'; '0136';...  
            '0137'; '0138'; '0140'; '0158'; '0162'; '0163'; '0178';...
            '0179'; '0255'; '0257'; '0348'; '0378'; '0384']; 
        
SUBJ_ASD = ['0106'; '0107'; '0139'; '0141'; '0159'; '0160'; '0161';...  
            '0164'; '0253'; '0254'; '0256'; '0273'; '0274'; '0275';...
            '0276'; '0346'; '0347'; '0351'; '0358'; '0357';...
            '0380'; '0381'; '0382'; '0383'];  
%without  0349; 0135; '0135'; 0379; 
SUBJ = [SUBJ_ASD; SUBJ_NT];

%loop for all subjects
for s=1: size (SUBJ_NT,1)
    subj = SUBJ_NT (s,:);
    
    load(strcat('/net/server/data/Archive/aut_gamma/orekhova/KI/Results_Alpha_and_Gamma/', subj, '/DICS_6conditions/', subj, '_alpha_sensor_spectr.mat'));
    savemegto = strcat(savepath, subj); 
    
    % calculate individual alpha freq for interstimulus mag
    alpha_slow_mag(s) = alpha_freq(SPECTR.freq,SPECTR.isi_mag{1}); 
    alpha_medium_mag(s) = alpha_freq(SPECTR.freq,SPECTR.isi_mag{2});
    alpha_fast_mag(s) = alpha_freq(SPECTR.freq,SPECTR.isi_mag{3});
    
    % calculate individual alpha freq for interstimulus grad 
    alpha_slow_grad(s) = alpha_freq(SPECTR.freq,SPECTR.isi_grad{1}); 
    alpha_medium_grad(s) = alpha_freq(SPECTR.freq,SPECTR.isi_grad{2}); 
    alpha_fast_grad(s) = alpha_freq(SPECTR.freq,SPECTR.isi_grad{3});
    
    % calculate individual alpha freq for stimulation mag
    alpha_slow_mag_post(s) = alpha_freq(SPECTR.freq,SPECTR.post_mag{1}); 
    alpha_medium_mag_post(s) = alpha_freq(SPECTR.freq,SPECTR.post_mag{2});
    alpha_fast_mag_post(s) = alpha_freq(SPECTR.freq,SPECTR.post_mag{3});
    
    % calculate individual alpha power for stimulation grad
    alpha_slow_grad_post(s) = alpha_freq(SPECTR.freq,SPECTR.post_grad{1}); 
    alpha_medium_grad_post(s) = alpha_freq(SPECTR.freq,SPECTR.post_grad{2});
    alpha_fast_grad_post(s) = alpha_freq(SPECTR.freq,SPECTR.post_grad{3});
    
    
    power_slow_mag(s) = alpha_power(SPECTR.freq, SPECTR.isi_mag{1}); 
    power_medium_mag(s) = alpha_power(SPECTR.freq, SPECTR.isi_mag{2}); 
    power_fast_mag(s) = alpha_power(SPECTR.freq, SPECTR.isi_mag{3}); 
    
    power_slow_grad(s) = alpha_power(SPECTR.freq, SPECTR.isi_grad{1}); 
    power_medium_grad(s) = alpha_power(SPECTR.freq, SPECTR.isi_grad{2}); 
    power_fast_grad(s) = alpha_power(SPECTR.freq, SPECTR.isi_grad{3}); 
    
    power_slow_mag_post(s) = alpha_power(SPECTR.freq, SPECTR.post_mag{1}); 
    power_medium_mag_post(s) = alpha_power(SPECTR.freq, SPECTR.post_mag{2}); 
    power_fast_mag_post(s) = alpha_power(SPECTR.freq, SPECTR.post_mag{3}); 
    
    power_slow_grad_post(s) = alpha_power(SPECTR.freq, SPECTR.post_grad{1}); 
    power_medium_grad_post(s) = alpha_power(SPECTR.freq, SPECTR.post_grad{2}); 
    power_fast_grad_post(s) = alpha_power(SPECTR.freq, SPECTR.post_grad{3}); 
    
end



% combine results for three conditions
isi_mag = [alpha_slow_mag', alpha_medium_mag', alpha_fast_mag', power_slow_mag', power_medium_mag', power_fast_mag'];
isi_grad = [alpha_slow_grad', alpha_medium_grad', alpha_fast_grad', power_slow_grad', power_medium_grad', power_fast_grad'];
post_mag = [alpha_slow_mag_post', alpha_medium_mag_post', alpha_fast_mag_post', power_slow_mag_post', power_medium_mag_post', power_fast_mag_post'];
post_grad = [alpha_slow_grad_post', alpha_medium_grad_post', alpha_fast_grad_post', power_slow_grad_post', power_medium_grad_post', power_fast_grad_post'];


%save as csv
filename = strcat(savepath, '1_results/individual_alpha_isi_mag_NT.csv');
csvwrite(filename, isi_mag);
filename = strcat(savepath, '1_results/individual_alpha_isi_grad_NT.csv');
csvwrite(filename, isi_grad);
filename = strcat(savepath, '1_results/individual_alpha_post_mag_NT.csv');
csvwrite(filename, post_mag);
filename = strcat(savepath, '1_results/individual_alpha_post_grad_NT.csv');
csvwrite(filename, post_grad);