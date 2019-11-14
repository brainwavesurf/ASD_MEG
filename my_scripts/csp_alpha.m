%% CSP analysis by fieldtrip function ft_componentanalysis

clear all;
close all;

fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(fieldtripfolder, path);
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);

realdatapath = '/home/a_shishkina/data/KI/SUBJECTS/';
savepath = '/home/a_shishkina/data/KI/Results_Alpha_and_Gamma/';

SUBJ_NT = [ '0101'; '0102'; '0103'; '0104'; '0105'; '0135'; '0136';...  
            '0137'; '0138'; '0140'; '0158'; '0162'; '0163'; '0178';...
            '0179'; '0255'; '0257'; '0348'; '0378'; '0379'; '0384']; 
        
SUBJ_ASD = ['0106'; '0107'; '0139'; '0141'; '0159'; '0160'; '0161';...  
            '0164'; '0253'; '0254'; '0256'; '0273'; '0274'; '0275';...
            '0276'; '0346'; '0347'; '0349'; '0351'; '0357'; '0358';...
            '0380'; '0381'; '0382'; '0383'];  
        
SUBJ = [SUBJ_NT; SUBJ_ASD];

SUBJ = ['0106'];

for s=1: size (SUBJ,1)
    close all
    subj = SUBJ (s,:); 
    savemegto = strcat(savepath, subj);
    epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
    
    %load alpha epochs
    load(strcat(epofolder, subj, '_preproc_alpha_epochs_10_13.mat'));
    
    cfg = [];
    data = ft_appenddata(cfg, slow_alpha_epochs, fast_alpha_epochs); %append two structural data

%%example for subj 0106
%
%     label: {306×1 cell}
%      trialinfo: [116×1 double]
%     sampleinfo: [116×2 double]
%           grad: [1×1 struct]
%          trial: {1×116 cell}
%           time: {1×116 cell}
%        fsample: 500
%            cfg: [1×1 struct]
    
    % prepare vectors that assigns slow and fast trials to class 1 or 2
    slow_label = zeros(size(slow_alpha_epochs.trialinfo)); slow_label(:) = 1;
    fast_label = zeros(size(fast_alpha_epochs.trialinfo)); fast_label(:) = 2;
    
    
    % The csp method implements the common-spatial patterns method
    cfg = [];
    cfg.channel = 'MEGMAG';
    cfg.method = 'csp';
    cfg.csp.classlabels = [slow_label; fast_label]; % vector that assigns a trial to class 1 or 2
    cfg.csp.numfilters  = 102; % the number of spatial filters to use
    [comp] = ft_componentanalysis(cfg, data);

%%example for subj 0106
%
%        fsample: 500
%           time: {1×116 cell}
%          trial: {1×116 cell}
%           topo: [102×87 double]
%       unmixing: [87×102 double]
%          label: {10×1 cell}
%      topolabel: {102×1 cell}
%           grad: [1×1 struct]
%     sampleinfo: [116×2 double]
%      trialinfo: [116×1 double]
%            cfg: [1×1 struct]

    filename = strcat(epofolder, subj, '_CSP_alpha_mag.mat');
    save(filename, 'comp');
end


figure(1)
cfg = [];
cfg.component = 1:30;       % the component(s) that should be plotted
cfg.layout    = 'neuromag306mag.lay'; % the layout file that should be used for plotting
cfg.comment   = 'no';
ft_topoplotIC(cfg, comp)

figure(2)
cfg = [];
cfg.component = 31:60;       % the component(s) that should be plotted
cfg.layout    = 'neuromag306mag.lay'; % the layout file that should be used for plotting
cfg.comment   = 'no';
ft_topoplotIC(cfg, comp)

figure(3)
cfg = [];
cfg.component = 61:length(comp.label);       % the component(s) that should be plotted
cfg.layout    = 'neuromag306mag.lay'; % the layout file that should be used for plotting
cfg.comment   = 'no';
ft_topoplotIC(cfg, comp)

figure(4)
cfg = [];
cfg.component = 1:6;       % the component(s) that should be plotted
cfg.layout    = 'neuromag306mag.lay'; % the layout file that should be used for plotting
cfg.comment   = 'no';
ft_topoplotIC(cfg, comp)

saveas(figure(1),[savepath, '/1_results/', '0106_ASD_1_30_components_mag.jpeg']);
saveas(figure(2),[savepath, '/1_results/', '0106_ASD_31_60_components_mag.jpeg']);
saveas(figure(3),[savepath, '/1_results/', '0106_ASD_61_89_components_mag.jpeg']);
saveas(figure(4),[savepath, '/1_results/', '0106_ASD_6_components_mag.jpeg']);

