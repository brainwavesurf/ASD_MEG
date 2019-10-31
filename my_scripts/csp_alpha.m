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

for s=1: size (SUBJ,1)
    close all
    subj = SUBJ (s,:); 
    savemegto = strcat(savepath, subj);
    epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
    
    load(strcat(epofolder, subj, '_preproc_alpha_epochs.mat'));
    cfg = [];
    data = ft_appenddata(cfg, slow_alpha_epochs, fast_alpha_epochs); %append two structural data
    
%     label: {306×1 cell}
%      trialinfo: [116×1 double]
%     sampleinfo: [116×2 double]
%           grad: [1×1 struct]
%          trial: {1×116 cell}
%           time: {1×116 cell}
%        fsample: 500
%            cfg: [1×1 struct]
    
    % prepare vectors that assigns slow and fast trials to class 1 or 2
    slow_label = zeros(59, 1); slow_label(:) = 1;
    fast_label = zeros(57, 1); fast_label(:) = 2;
    
    % The csp method implements the common-spatial patterns method
    cfg = [];
    cfg.method = 'csp';
    cfg.csp.classlabels = [slow_label; fast_label]; % vector that assigns a trial to class 1 or 2
    cfg.csp.numfilters  = 10; % the number of spatial filters to use
    [comp] = ft_componentanalysis(cfg, data);
    
%        fsample: 500
%           time: {1×116 cell}
%          trial: {1×116 cell}
%           topo: [306×10 double]
%       unmixing: [10×306 double]
%          label: {10×1 cell}
%      topolabel: {306×1 cell}
%           grad: [1×1 struct]
%     sampleinfo: [116×2 double]
%      trialinfo: [116×1 double]
%            cfg: [1×1 struct]
end

figure(1)
cfg = [];
cfg.component = 1:10;       % the component(s) that should be plotted
cfg.layout    = 'neuromag306mag.lay'; % the layout file that should be used for plotting
cfg.comment   = 'no';
ft_topoplotIC(cfg, comp)
saveas(figure(1),[savepath, '/1_results/', 'CSP_components.jpeg']);
