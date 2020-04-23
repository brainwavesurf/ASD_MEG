% comments for DATA from ONE SUBJ #0107 in 2 conditions and in an epoch of 401 sampling points, srate = 500Hz, MEGMAG data
% ntrial_slow = 64 trials, ntrial_fast = 65
close all;
tStart = tic;
fieldtripfolder = '/home/a_shishkina/fieldtrip/';
path(fieldtripfolder, path);
ft_defaults;
path('/home/a_shishkina/fieldtrip/external/mne/', path);
realdatapath = '/net/server/data/Archive/aut_gamma/orekhova/KI/SUBJECTS/';
savepath = '/net/server/data/Archive/aut_gamma/orekhova/KI/Scripts_bkp/Shishkina/KI/Results_Alpha_and_Gamma/';


SUBJ_NT = [ '0101'; '0102'; '0103'; '0104'; '0105'; '0135'; '0136';...  
            '0137'; '0138'; '0140'; '0158'; '0162'; '0163'; '0178';...
            '0179'; '0255'; '0257'; '0348'; '0378'; '0379'; '0384']; 
        
SUBJ_ASD = ['0106'; '0107'; '0139'; '0141'; '0159'; '0160'; '0161';...  
            '0164'; '0253'; '0254'; '0256'; '0273'; '0274'; '0275';...
            '0276'; '0346'; '0347'; '0349'; '0351'; '0358';...
            '0380'; '0381'; '0382'; '0383'];  

SUBJ = [SUBJ_NT; SUBJ_ASD];

for s=1: size (SUBJ,1)
    close all
    subj = SUBJ (s,:); 
    
    %load alpha epochs
    load(strcat(savepath, subj, '/', subj, '_preproc_alpha_10_17_epochs.mat'));

    cfg = [];
    cfg.channel = 'MEGMAG';
    data_slow = ft_selectdata(cfg, slow_alpha_isi); 
    data_fast = ft_selectdata(cfg, fast_alpha_isi); 

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

%     figure;plot(diag(D1), '.-k')
%     xlabel('#eigenvalue')
%     ylabel('eigenvalue')
%     saveas(figure(2), [savepath, '1_results/CSP_matlab_plots/', subj, '_eigenvalue.jpeg']); 

    A1 = inv(W1); % nIC x nMAG: 102x102; filters must be positive and negative: ok 

    indsel = 1:3; %pick 3 eigenvalues 
    W1 = W1(:,indsel); %nMAG x nCSP: 102x3 , mixing matrix
    A1 = A1(indsel,:); %nCSP x nMAG: 3x102 , unmixing matrix

    Xcsp_slow = nan(ntrial_slow, sampoints_slow, size(W1,2)); %64x401x3
    Xcsp_fast = nan(ntrial_fast, sampoints_fast, size(W1,2)); %65x401x3
    
    %% pattern maximizing the differences between cov1/cov2

    pattern_ICcsp_slowVSfast = transpose(A1); % A_CSP [nCSP x nMEG] 102 x 3

    for j = 1:ntrial_slow   
    Xcsp_slow(j,:,:)=transpose(squeeze(data_slow_tot(:,:,j)))*W1; % ntrial x nsampl x 3 = [ntrial x nsampl x nIC] * [nIC x 3] 
    end  

    for j = 1:ntrial_fast  
    Xcsp_fast(j,:,:)=transpose(squeeze(data_fast_tot(:,:,j)))*W1; %% 3 x nsampl x ntrial = [3 x nIC ] x [nIC x nsampl x ntrial]  
    end

    filename = strcat(savepath, subj, '/', subj, '_csp_analysis_1_3.mat');
    save(filename, 'W1', 'A1', 'pattern_ICcsp_slowVSfast', 'Xcsp_fast', 'Xcsp_slow');
    
end 
tEnd = toc(tStart);