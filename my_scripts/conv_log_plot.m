%plot log power after mtmconvol for posterior gradiometers for all subjects in fast and slow conditions

clear all;
close all;

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(path, fieldtripfolder)
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);

savepath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';
realdatapath = '/home/a_shishkina/data/KI/SUBJECTS/';

%load subj list
SUBJ = ['0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0136'; '0137'; '0138'; '0139'; '0140'; '0141'; '0158'; '0159'; '0160'; '0161'; '0162'; '0163'; '0164'; '0178'; '0179'; '0253'; '0254'; '0255'; '0256'; '0257'; '0259'; '0273'; '0274'; '0275'; '0276'; '0277'; '0346'; '0347'; '0348'; '0350'; '0351'; '0357'; '0358'; '0378'; '0380'; '0381'; '0382'; '0383'; '0384']; 
%load posterior sensors list
post_sens = {'MEG1932',  'MEG1922', 'MEG2042',  'MEG2032',  'MEG2112', 'MEG2122',  'MEG2342', 'MEG2332',  'MEG1732', 'MEG1942', 'MEG1912', 'MEG2012', 'MEG2022', 'MEG2312', 'MEG2322', 'MEG2512',...
             'MEG1933',  'MEG1923', 'MEG2043',  'MEG2033',  'MEG2113', 'MEG2123',  'MEG2343', 'MEG2333',  'MEG1733', 'MEG1943', 'MEG1913', 'MEG2013', 'MEG2023', 'MEG2313', 'MEG2323', 'MEG2513'};

%loop for all subjects
for s=1: size (SUBJ,1)
    
    subj = SUBJ (s,:); 
    savemegto = strcat(savepath, subj);
    %load freqanalysis results
    freq_all = load(strcat(savepath, subj, '/', subj, '_freqanalysis.mat'));
    %select mtmconv for fast and slow conditions
    freqFast_grad  = freq_all.wvlts_fast_grad;
    freqSlow_grad  = freq_all.wvlts_slow_grad;
    %select posterior sensors
    cfg = [];
    cfg.channel = post_sens;
    freqFast_postsens{s}  = ft_selectdata(cfg, freqFast_grad);
    freqSlow_postsens{s}  = ft_selectdata(cfg, freqSlow_grad);
    
%        label: {32×1 cell}
%          freq: [1×26 double]
%          time: [1×45 double]
%     powspctrm: [32×26×45 double]
%          grad: [1×1 struct]
%           cfg: [1×1 struct]
%        dimord: 'chan_freq_time'
end


%do averaging over channels 
cfg = [];
cfg.avgoverchan = 'yes';
fast_avg_grad = ft_selectdata(cfg, freqFast_postsens{:});
slow_avg_grad = ft_selectdata(cfg, freqSlow_postsens{:});

%save stats
filename = strcat(savepath, '1_results/', 'fast_slow_avg_grad.mat');
save(filename, 'fast_avg_grad', 'slow_avg_grad');

%  label: {'mean(MEG1732, MEG1733, MEG1912, MEG1913, MEG1922, MEG1923, MEG1932, MEG1933, MEG1942, MEG1943, MEG2012, MEG2013, MEG2022, MEG2023, MEG2032, MEG2033, MEG2042, MEG2043, MEG2112, MEG2113, MEG2122, MEG2123, MEG2312, MEG2313, MEG2322, MEG2323, MEG2332, MEG2333, MEG2342, MEG2343, MEG2512, MEG2513)'}
%          freq: [1×26 double]
%          time: [1×45 double]
%     powspctrm: [1×26×45 double]
%          grad: [1×1 struct]
%           cfg: [1×1 struct]
%        dimord: 'chan_freq_time'


%plot figures
figure(1)
subplot(2,1,1)
plot(slow_avg_grad.freq, slow_avg_grad.powspctrm(1,:,11), '-g'); title('FFT power, avg posterior sens for all subj')
hold on
plot(fast_avg_grad.freq, fast_avg_grad.powspctrm(1,:,11), '-b'); 
legend('slow', 'fast'); xlim([5 30]);

subplot(2,1,2)
plot(slow_avg_grad.freq, log(slow_avg_grad.powspctrm(1,:, 11)), '-g'); title('FFT log power, avg posterior sens for all subj')
hold on
plot(fast_avg_grad.freq, log(fast_avg_grad.powspctrm(1,:, 11)), '-b');
legend('slow', 'fast'); xlim([5 30]);

saveas(figure(1),[savepath, '/1_results/', 'conv_plot_prestim.jpeg']);
