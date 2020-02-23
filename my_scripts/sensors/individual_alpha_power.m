%calculation individual alpha freq with linear trend removal
clear all;
close all;

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(path, fieldtripfolder)
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);
path('/home/a_shishkina/projects/ASD_MEG/Git/ASD_MEG/my_scripts', path);
path('/home/a_shishkina/projects/ASD_MEG/Git/ASD_MEG/my_scripts/sensors',path);
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

%matrix for function output
output_freq = zeros(size(SUBJ,1),12,4);
output_pow = zeros(size(SUBJ,1),12,2);
output_pow_cluster = zeros(size(SUBJ,1),12,2);

%indexes for different conditions and different channels
isi_mag_idx = [1:3];
isi_grad_idx = [4:6];
post_mag_idx = [7:9];
post_grad_idx = [10:12];

%indexes for different values
with_trend_idx = 1;
without_trend_idx = 2;
slope_idx = 3;
intercept_idx = 4;

%loop for all subjects
for s=1:size(SUBJ,1) 
    
    subj = SUBJ(s,:); 
    load(strcat('/net/server/data/Archive/aut_gamma/orekhova/KI/Results_Alpha_and_Gamma/', subj, '/DICS_6conditions/', subj, '_alpha_sensor_spectr.mat'));
    savemegto = strcat(savepath, subj); 
    
    x = SPECTR.freq;
    y = [SPECTR.isi_mag{1};SPECTR.isi_mag{2};SPECTR.isi_mag{3};...
        SPECTR.isi_grad{1};SPECTR.isi_grad{2};SPECTR.isi_grad{3};...
        SPECTR.post_mag{1};SPECTR.post_mag{2};SPECTR.post_mag{3};...
        SPECTR.post_grad{1};SPECTR.post_grad{2};SPECTR.post_grad{3}];

    for i = 1:12 %through y 
        [output_freq(s,i,with_trend_idx), output_freq(s,i,without_trend_idx), output_freq(s,i,slope_idx), output_freq(s,i,intercept_idx)] = alpha_freq(x, y(i,:));
        [output_pow(s,i,with_trend_idx), output_pow(s,i,without_trend_idx)] = alpha_power(x, y(i,:));
        [output_pow_cluster(s,i,with_trend_idx), output_pow_cluster(s,i,without_trend_idx)] = alpha_power_cluster(x, y(i,:));
    end
end

%combine results 
isi_mag     = [output_freq(:,isi_mag_idx,with_trend_idx), output_freq(:,isi_mag_idx,without_trend_idx), output_freq(:,isi_mag_idx,slope_idx), output_freq(:,isi_mag_idx,intercept_idx), ...
               output_pow(:,isi_mag_idx,with_trend_idx), output_pow(:,isi_mag_idx,without_trend_idx), output_pow_cluster(:,isi_mag_idx,with_trend_idx), output_pow_cluster(:,isi_mag_idx,without_trend_idx)];

isi_grad    = [output_freq(:,isi_grad_idx,with_trend_idx), output_freq(:,isi_grad_idx,without_trend_idx), output_freq(:,isi_grad_idx,slope_idx), output_freq(:,isi_grad_idx,intercept_idx), ...
               output_pow(:,isi_grad_idx,with_trend_idx), output_pow(:,isi_grad_idx,without_trend_idx), output_pow_cluster(:,isi_grad_idx,with_trend_idx), output_pow_cluster(:,isi_grad_idx,without_trend_idx)];

post_mag    = [output_freq(:,post_mag_idx,with_trend_idx), output_freq(:,post_mag_idx,without_trend_idx), output_freq(:,post_mag_idx,slope_idx), output_freq(:,post_mag_idx,intercept_idx), ...
               output_pow(:,post_mag_idx,with_trend_idx), output_pow(:,post_mag_idx,without_trend_idx), output_pow_cluster(:,post_mag_idx,with_trend_idx), output_pow_cluster(:,post_mag_idx,without_trend_idx)];

post_grad   = [output_freq(:,post_grad_idx,with_trend_idx), output_freq(:,post_grad_idx,without_trend_idx), output_freq(:,post_grad_idx,slope_idx), output_freq(:,post_grad_idx,intercept_idx), ...
               output_pow(:,post_grad_idx,with_trend_idx), output_pow(:,post_grad_idx,without_trend_idx), output_pow_cluster(:,post_grad_idx,with_trend_idx), output_pow_cluster(:,post_grad_idx,without_trend_idx)];

%save as csv 
filename = strcat(savepath, '1_results/individual_alpha_isi_mag.csv');
csvwrite(filename, isi_mag);
filename = strcat(savepath, '1_results/individual_alpha_isi_grad.csv');
csvwrite(filename, isi_grad);
filename = strcat(savepath, '1_results/individual_alpha_post_mag.csv');
csvwrite(filename, post_mag);
filename = strcat(savepath, '1_results/individual_alpha_post_grad.csv');
csvwrite(filename, post_grad);
