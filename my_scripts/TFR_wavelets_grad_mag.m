%TFR with wavelets for mag and grad
%plot TFR

SUBJ = ['0076'; '0101'; '0102'; '0103'; '0104'; '0105'; '0106'; '0107'; '0136'; '0137'; '0138'; '0139'; '0140'; '0141'; '0158'; '0159'; '0160'; '0161'; '0162'; '0163'; '0164'; '0178'; '0179'; '0253'; '0254'; '0255'; '0256'; '0257'; '0259'; '0273'; '0274'; '0275'; '0276'; '0277'; '0346'; '0347'; '0348'; '0350'; '0351'; '0357'; '0358'; '0378'; '0380'; '0381'; '0382'; '0383'; '0384']; 

for s=1: size (SUBJ,1)
    
    subj = SUBJ (s,:); 
    savemegto = strcat(savepath, subj);
    epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
    
    epo = load(strcat(epofolder, subj, '_preproc_epochs.mat'));
    cfg = [];
    cfg.channel = 'MEGMAG';
    epo_fast_mag = ft_selectdata(cfg, epo.fast_epochs);
    epo_slow_mag = ft_selectdata(cfg, epo.slow_epochs);
    
    cfg.channel = 'MEGGRAD';
    epo_fast_grad = ft_selectdata(cfg, epo.fast_epochs);
    epo_slow_grad = ft_selectdata(cfg, epo.slow_epochs);

%% Calculate the TFRs for the two experimental conditions (Fast and Slow) for each subject
    cfg = [];
    cfg.channel    = 'MEG';
    cfg.method     = 'wavelet';
    cfg.width      = 3;
    cfg.output     = 'pow';
    cfg.foi        = 5:20;
    cfg.toi        = -1:0.05:1.2;

    freqFast_mag{s}  = ft_freqanalysis(cfg, epo_fast_mag);
    freqSlow_mag{s}  = ft_freqanalysis(cfg, epo_slow_mag);
    
    freqFast_grad{s}  = ft_freqanalysis(cfg, epo_fast_grad);
    freqSlow_grad{s}  = ft_freqanalysis(cfg, epo_slow_grad);
    
end
%% Calculate the grand averages of the TFRs for the fast and slow conditions

cfg = [];
cfg.foilim = [5 20];
cfg.toilim = [-1 1.2];
slow_avg_mag = ft_freqgrandaverage(cfg, freqSlow_mag{:});
fast_avg_mag = ft_freqgrandaverage(cfg, freqFast_mag{:});

slow_avg_grad = ft_freqgrandaverage(cfg, freqSlow_grad{:});
fast_avg_grad = ft_freqgrandaverage(cfg, freqFast_grad{:});

%Averaging of grads info because of warning "discarding gradiometer information because it cannot be averaged" 
input_fast_mag = []; input_slow_mag = [];
input_fast_grad = []; input_slow_grad = [];

for i = 1:size (SUBJ,1)
    input_fast_mag  = [input_fast_mag, freqFast_mag{i}.grad];
    input_slow_mag  = [input_slow_mag, freqSlow_mag{i}.grad];
   
end

input_fast_grad = [input_slow_grad, freqFast_grad{i}.grad];
input_slow_grad = [input_slow_grad, freqSlow_grad{i}.grad];

AvgGrad_fast_mag = ft_average_sens(input_fast_mag);
fast_avg_mag.grad = AvgGrad_fast_mag;

AvgGrad_slow_mag = ft_average_sens(input_slow_mag);
slow_avg_mag.grad = AvgGrad_slow_mag;

AvgGrad_fast_grad = ft_average_sens(input_fast_grad);
fast_avg_grad.grad = AvgGrad_fast_grad;

AvgGrad_slow_grad = ft_average_sens(input_slow_grad);
slow_avg_grad.grad = AvgGrad_slow_grad;

% Plot TFR wavelet for averaged data
figure(1);
cfg = [];
cfg.baseline     = [0 inf];
cfg.baselinetype = 'absolute';
cfg.maskstyle    = 'saturation';
cfg.zlim         = 'absmax';
cfg.channel      = 'MEG2031';
cfg.layout       = 'neuromag306mag.lay';
cfg.marker       = 'on';
subplot(2,1,1); ft_singleplotTFR(cfg, slow_avg_mag); title('width = 3, slow, mag'); 
subplot(2,1,2); ft_singleplotTFR(cfg, fast_avg_mag); title('width = 3, fast, mag');
saveas(figure(1), [savepath, '1_results/', 'avg_TFwavelet_width5_mag.jpeg']);

figure(2);
cfg.channel      = 'MEG2043';
cfg.layout       = 'neuromag306planar.lay';
subplot(2,1,1); ft_singleplotTFR(cfg, slow_avg_grad); title('width = 3, slow, grad'); 
subplot(2,1,2); ft_singleplotTFR(cfg, fast_avg_grad); title('width = 3, fast, grad');
saveas(figure(2), [savepath, '1_results/', 'avg_TFwavelet_width3_grad.jpeg']);
