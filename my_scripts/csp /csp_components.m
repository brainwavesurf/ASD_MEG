% comments for DATA from ONE SUBJ #0107 in 2 conditions and in an epoch of 401 sampling points, srate = 500Hz, MEGMAG data
% ntrial_slow = 64 trials, ntrial_fast = 65           

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
            '0276'; '0346'; '0347'; '0349'; '0351'; '0358';...
            '0380'; '0381'; '0382'; '0383'];  
%without '0357';
SUBJ = [SUBJ_NT; SUBJ_ASD];

for s=1: size (SUBJ,1)
   
    subj = SUBJ (s,:); 
    savemegto = strcat(savepath, subj);
    epofolder = strcat(realdatapath, subj, '/ICA_nonotch_crop', '/epochs/');
    
    %load alpha epochs
    load(strcat(epofolder, subj, '_preproc_alpha_bp_epochs.mat'));

    cfg = [];
    cfg.channel = 'MEGMAG';
    data_slow = ft_selectdata(cfg, slow_alpha_bp); %slow_alpha_post
    data_fast = ft_selectdata(cfg, fast_alpha_bp); %fast_alpha_post

    ntrial_slow = size(data_slow.trial,2); %number of trials in the slow condition: 64
    ntrial_fast = size(data_fast.trial,2); %number of trials in the fast condition: 65

    sampoints_slow = size(data_slow.time{1},2); %number of sampling point in slow and fast conditions: 401
    sampoints_fast = size(data_fast.time{1},2);

    data_slow_tot = nan(102,sampoints_slow,ntrial_slow); % nMAG x nsamp x ntrial: 102x401x64 
    data_fast_tot = nan(102,sampoints_fast,ntrial_fast); % nMAG x nsamp x ntrial: 102x401x65 
    
    cov_slow = zeros(102,102,ntrial_slow); % nMAG x nIC x ntrial: 102x102x64
    cov_fast = zeros(102,102,ntrial_fast); % nMAG x nIC x ntrial: 102x102x65

    %concatenate epochs
    for k = 1 : ntrial_slow
    data_slow_tot(:,:,k) = data_slow.trial{k};
    cov_slow(:,:,k)= cov(transpose(squeeze(data_slow_tot(:,:,k)))); 
    end

    for k = 1 : ntrial_fast
    data_fast_tot(:,:,k) = data_fast.trial{k};
    cov_fast(:,:,k)= cov(transpose(squeeze(data_fast_tot(:,:,k))));
    end

    cov_slow = mean(cov_slow,3); %102x102
    cov_fast = mean(cov_fast,3); %102x102         
    %% select 3 eigenvectors from each tail of eigenvalues [W1, D1] = eig(cov_slow_cond, cov_slow_cond + cov_fast_cond);
    % eig(cov1, cov2) is similar, but different normalization

    %           Condition1  Condition1+Condition2
    [W1,D1] = eig(cov_slow, cov_slow + cov_fast); %W1: spatial filters

    figure;plot(diag(D1),'m*')
    
    A1 = inv(W1); % nIC x nMAG: 102x102; filters must be positive and negative: ok 

    indsel = [1:3, size(W1,1)-2:size(W1,1)]; %pick 3 eigenvalues from each tail
    W1 = W1(:,indsel); %nMAG x nCSP: 102x6 , mixing matrix
    A1 = A1(indsel,:); %nCSP x nMAG: 6x102 , unmixing matrix

    Xcsp_slow = nan(ntrial_slow, sampoints_slow, size(W1,2)); %64x401x6
    Xcsp_fast = nan(ntrial_fast, sampoints_fast, size(W1,2)); %65x401x6
    
    %% pattern maximizing the differences between cov1/cov2

    pattern_ICcsp_slowVSfast = transpose(A1); % A_CSP [nCSP x nMEG] 102 x 6

    for j = 1:ntrial_slow   
    Xcsp_slow(j,:,:)=transpose(squeeze(data_slow_tot(:,:,j)))*W1; % ntrial x nsampl x 6 = [ntrial x nsampl x nIC] * [nIC x 6] 
    end  

    for j = 1:ntrial_fast  
    Xcsp_fast(j,:,:)=transpose(squeeze(data_fast_tot(:,:,j)))*W1; %% 6 x nsampl x ntrial = [6 x nIC ] x [nIC x nsampl x ntrial]  
    end

    filename = strcat(savepath, subj, '/', subj, '_csp_analysis.mat');
    save(filename, 'W1', 'A1', 'pattern_ICcsp_slowVSfast', 'Xcsp_fast', 'Xcsp_slow');
end   
    %% plot results 
    %oscillatory 
    for i = 1:6
    figure(1)
    subplot(2,3,i)
    plot(Xcsp_fast(1,:,i), 'r')
    hold on
    plot(Xcsp_slow(1,:,i), 'b')
    legend('fast','slow')
    xlim([0 401])
    title(['compoment', num2str(i)])
    end
    
    %topo
    slow_label = zeros(ntrial_slow, 1); slow_label(:) = 1;
    fast_label = zeros(ntrial_fast, 1); fast_label(:) = 2;

    data = ft_appenddata(cfg, data_slow, data_fast); %append two structural data

    cfg = [];
    cfg.method = 'csp';
    cfg.csp.classlabels = [slow_label; fast_label]; % vector that assigns a trial to class 1 or 2
    cfg.csp.numfilters  = 10; % the number of spatial filters to use
    [comp] = ft_componentanalysis(cfg, data);
    comp.topo = pattern_ICcsp_slowVSfast;
    comp.mixing = W1;
    comp.unmixing = A1;

    figure(2)
    cfg = [];
    cfg.component = 1:6;       % the component(s) that should be plotted
    cfg.layout    = 'neuromag306mag.lay'; % the layout file that should be used for plotting
    cfg.comment   = 'no';
    ft_topoplotIC(cfg, comp)

%     filename = strcat(savepath, subj, '/', subj, '_csp_components.jpeg');
%     saveas(figure(1), filename);
%     filename = strcat(savepath, subj, '/', subj, '_csp_components_topo.jpeg');
%     saveas(figure(2), filename);
end
% 
