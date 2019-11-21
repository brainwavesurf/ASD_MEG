% DATA from ONE SUBJ #0106 in 2 conditions and in an epoch of 401 sampling points, srate = 500Hz, MEGMAG data
% ntrial_slow = 59 trials, ntrialslow = 57           

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
    cfg.channel = 'MEGMAG';
    data_slow = ft_selectdata(cfg, slow_alpha_epochs);
    data_fast = ft_selectdata(cfg, fast_alpha_epochs);

    ntrial_slow = size(slow_alpha_epochs.trial,2); %number of trials in the slow condition: 59
    ntrial_fast = size(fast_alpha_epochs.trial,2); %number of trials in the fast condition: 57

    sampoints_slow = size(data_slow.time{1},2); %number of sampling point in slow and fast conditions: 401
    sampoints_fast = size(data_fast.time{1},2);

    data_slow_tot = nan(102,sampoints_slow,ntrial_slow); % nMAG x nsamp x ntrial: 102x401x59 
    data_fast_tot = nan(102,sampoints_fast,ntrial_fast); % nMAG x nsamp x ntrial: 102x401x57 

    cov_slow = zeros(102,102,ntrial_slow); % nMAG x nIC x ntrial: 102x102x59
    cov_fast = zeros(102,102,ntrial_fast); % nMAG x nIC x ntrial: 102x102x57

    %concatenate epochs
    for k = 1 : ntrial_slow
    data_slow_tot(:,:,k) = data_slow.trial{k};
    cov_slow(:,:,k)= cov(transpose(squeeze(data_slow_tot(:,:,k)))); 
    end

    for k = 1 : ntrial_fast
    data_fast_tot(:,:,k) = data_fast.trial{k};
    cov_fast(:,:,k)= cov(transpose(squeeze(data_fast_tot(:,:,k))));
    end

    cov_slow_cond = mean(cov_slow,3); %102x102
    cov_fast_cond = mean(cov_fast,3); %102x102         
    %% select 3 eigenvectors from each tail of eigenvalues [W1, D1] = eig(cov_slow_cond, cov_slow_cond + cov_fast_cond);

    [W1,D1] = eig(cov_slow_cond, cov_slow_cond + cov_fast_cond); %W1: spatial filters

    figure;plot(diag(D1),'m*')
    saveas(figure(1),[savepath, '/1_results/', '0106_ASD_eigenvalues_mag.jpeg'])

    mixing = W1; % nMAG x nIC : 102x102
    A1 = inv(W1); % nIC x nMAG: 102x102; filters must be positive and negative: ok 

    indsel = [1:3, size(W1,1)-2:size(W1,1)]; %pick 3 eigenvalues from each tail
    W1 = W1(:,indsel); %nMAG x nCSP: 102x6 , mixing matrix
    A1 = A1(indsel,:); %nCSP x nMAG: 6x102 , unmixing matrix

    Xcsp_slow = nan(ntrial_slow, sampoints_slow, size(W1,2)); %59x401x6
    Xcsp_fast = nan(ntrial_fast, sampoints_fast, size(W1,2)); %57x401x6

    for j = 1:ntrial_slow   
    Xcsp_slow(j,:,:)=transpose(squeeze(data_slow_tot(:,:,j)))*W1; % ntrial x nsampl x 6 = [ntrial x nsampl x nIC] * [nIC x 6] 
    end  

    for j = 1:ntrial_fast  
    Xcsp_fast(j,:,:)=transpose(squeeze(data_fast_tot(:,:,j)))*W1; %% 6 x nsampl x ntrial = [6 x nIC ] x [nIC x nsampl x ntrial]  
    end

    %for each Xcsp, find MEG pattern
    Pattern = mixing * transpose(A1); %102 x 6: [nMAG x nIC] * [nMAG x nCSP]
    
    %% calculate csp by fieldtrip for plot
    %load alpha epochs
    load(strcat(epofolder, subj, '_preproc_alpha_epochs_10_13.mat'));
    
    cfg = [];
    data = ft_appenddata(cfg, slow_alpha_epochs, fast_alpha_epochs); %append two structural data
    
    % prepare vectors that assigns slow and fast trials to class 1 or 2
    slow_label = zeros(size(slow_alpha_epochs.trialinfo)); slow_label(:) = 1;
    fast_label = zeros(size(fast_alpha_epochs.trialinfo)); fast_label(:) = 2;
    
    
    % The csp method implements the common-spatial patterns method
    cfg = [];
    cfg.channel = 'MEGMAG';
    cfg.method = 'csp';
    cfg.csp.classlabels = [slow_label; fast_label]; % vector that assigns a trial to class 1 or 2
    cfg.csp.numfilters  = 6; % the number of spatial filters to use
    [comp] = ft_componentanalysis(cfg, data);
    
    filename = strcat(epofolder, subj, '_csp_ftcompanalysis.mat');
    save(filename, 'comp');
    
%     % replace topo by pattern and unmixing by A1
%     cfg = [];   
%     comp.topo = Pattern;
%     comp.unmixing = A1;
%     cfg.component = 1:6; % the component(s) that should be plotted
%     cfg.layout    = 'neuromag306mag.lay'; % the layout file that should be used for plotting
%     cfg.comment   = 'no';
%     ft_topoplotIC(cfg, comp)
% 
%     saveas(figure(1),[savepath, '/1_results/', '0106_ASD_6_csp_comp_mag_.jpeg']);
    
    %% time series after csp for 6 patterns 
    for j = 1:ntrial_slow   
        Xcsp_timeseries_slow(:,:,j) = Pattern*transpose(squeeze(Xcsp_slow(j,:,:))); % nMAG x nsampl x ntrial =  [nMAG x 6] * [6 x nsampl x ntrial]
        data_slow.trial{j} = Xcsp_timeseries_slow(:,:,j); 
    end  

    for j = 1:ntrial_fast   
        Xcsp_timeseries_fast(:,:,j) = Pattern*transpose(squeeze(Xcsp_fast(j,:,:))); % nMAG x nsampl x ntrial =  [nMAG x 6] * [6 x nsampl x ntrial]
        data_fast.trial{j} = Xcsp_timeseries_fast(:,:,j);
    end

    filename = strcat(epofolder, subj, '_csp_analysis.mat');
    save(filename, 'mixing', 'W1', 'A1', 'Pattern', 'Xcsp_timeseries_slow', 'Xcsp_timeseries_fast', 'Xcsp_fast', 'Xcsp_slow');
%     cfg = [];
%     cgf.channel = 'megmag';
%     cfg.layout    = 'neuromag306mag.lay'; % the layout file that should be used for plotting
%     figure(2)
%     ft_topoplotTFR(cfg, data_slow)
%     figure(3)
%     ft_topoplotTFR(cfg, data_fast)
%     
%     saveas(figure(2),[savepath, '/1_results/', '0106_csp_timeseries_mag_slow.jpeg']);
%     saveas(figure(3),[savepath, '/1_results/', '0106_csp_timeseries_mag_fast.jpeg']);
end