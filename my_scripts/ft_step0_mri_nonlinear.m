% Preparation of volumetric head model for beamformer

%%
clear;
close all;
clc;

gridres = 6; % grid step in mm
step = gridres/10; % grid step in cm

% NB: add to path: 
fieldtripfolder = '/home/kolai/Documents/Shishkina/ProgramFiles/Fieldtrip/';
path(path, fieldtripfolder)
ft_defaults;
path('/home/kolai/Documents/Shishkina/ProgramFiles/Fieldtrip/external/mne/', path);

% This is an external matlab package used to save figures to PPTX
%path(path,'/External_matlab/exportToPPTX-master/');


%% 

megpath   = '/home/kolai/Documents/Shishkina/NeuralDataAnalysis/Autism/0101/MEGdata/';
mripath   = '/home/kolai/Documents/Shishkina/NeuralDataAnalysis/Autism/0101/freesurfer/';
DATAPATH  = '/home/kolai/Documents/Shishkina/NeuralDataAnalysis/Autism/0101/MEGdata/FT_beamf/';  


templatefile = strcat (fieldtripfolder, 'external/spm8/templates/T1.nii');
template_mri = ft_read_mri(templatefile);

% %     % ft template grid:
% %     templatedir = '/Applications/fieldtrip-20170515/template/sourcemodel/'
% %     template = load(fullfile(templatedir, 'standard_sourcemodel3d6mm'));

% Segmentation can be done with different brain thresholds:
brainthreshold = 0.5;

% SUBJ = [ 'G001'; 'G001'; 'G003'; 'G003'; 'G005';'G006';'G007'; 'G007';  'G013'; 'G014'; 'G015'; 'G015'; 'G016'; 'G016'; 'G017'; 'G017'; 'G020'; 'G020'; 'G021'; 'G021'; 'G022'; 'G022'; 'G023'; 'G023'; 'G025'; 'G025'; 'G026'; 'G026'; 'G027'; 'G027'; 'G028'; 'G028'; 'G030'; 'G030'];
% DATE = [   'f';    'l';    'f';    'l';    'l';   'f';   'f';   'l'; 'l'; 'l'; 'f'; 'l'; 'f'; 'l'; 'f'; 'l';  'f'; 'l'; 'f'; 'l'; 'f'; 'l'; 'f'; 'l';'f'; 'l'; 'f'; 'l'; 'f'; 'l'; 'f'; 'l'; 'f'; 'l'];

SUBJ = ['0101'];
DATE = ['oct-6'];


%  
% 
%%
s=1;
close all
subj = '0101'; 
date = 'oct-6';
%mrifolder =strcat( '/mri_linwarp_', num2str(gridres), 'mm_', 'brthr', '0.5', '/');
mrifolder = strcat( '/mri_nonlinwarp_', num2str(gridres), 'mm_', 'brthr', '0.5', '/'); 
brainthreshold=0.5;

%%
    
subjfolder = strcat(DATAPATH, subj, '_', date);
mkdir (subjfolder);
savemrito = strcat(subjfolder,  mrifolder);
mkdir (savemrito);
    
%%
%start PPTX report
%exportToPPTX('new');
%% transpath = '/home/eorekhova/FT_beamf_PMS/trans/';
DBPath = '';
transfile = strcat(DATAPATH, subj, '_', date, '/', subj, '_', date,'trans.mat');  % coordinate transformation file saved from MNE and exported to matlab
mridata = strcat(mripath , 'Case',  subj, '/mri/T1.mgz');
%mridata ='/Users/mtw/MEG/Karolinska/MRI_nii_subj/0137/0137_T1_192_RAS.mgz';
dataset = strcat(megpath, subj, '_rings_ICA_raw.fif'); 
    
%% read mri
mri_orig = ft_read_mri(strcat(mridata)); %??? error "Compressed .mgz files cannot be read on a PC"
mri_orig = ft_convert_units(mri_orig, 'cm'); 
%mri_orig = ft_determine_coordsys(mri_orig, 'interactive', 'yes');

    
%% This was figured out by EO, this is not the standart way, but it works!
% -transformation correction. NB: take mri_orig.hdr.tkrvox2ras instead of 
% mri_orig.transform
% load -trans.mat, it is the MNE MEG-head/MRI-RAS transformation matrix
load(transfile) % MEG/MRI coordinate transformation
trans.trans(1:3,4) = trans.trans(1:3,4)*100;% translation: meters to cm
%mri_orig.transform = inv(trans.trans)*mri_orig.transform; %trans.trans is 
%MEG to MRI/RAS;  mri_orig.transform is ALS to RAS 
ttt = mri_orig.transform;    %mri_orig.hdr.tkrvox2ras; 
ttt= mri_orig.hdr.tkrvox2ras; % This is for FS T1.mgz!%
ttt(1:3,:)=ttt(1:3,:)/10;
mri_orig.transform = inv(trans.trans)*(ttt);
mri_orig = ft_determine_coordsys(mri_orig, 'interactive', 'no');
mri_orig.coordsys='neuromag';
%  T = ttt; inv(trans.trans)/ttt;
%  mri_orig = ft_transform_geometry(T, mri_orig);

%% plot MEG headshape and sensors
grad = ft_read_sens(dataset,'senstype','meg'); 
grad  = ft_convert_units(grad , 'cm');
ft_datatype_sens(grad)
shape   = ft_read_headshape(dataset);
shape = ft_convert_units(shape , 'cm');
h=figure;
ft_plot_headshape(shape)
ft_plot_sens(grad, 'style', '*b');
view([1 0 0])
title('MEG headshape and sensors', 'FontSize', 13)
%print -dpng natmeg_dip_geometry1.png
%exportToPPTX('addslide');
%exportToPPTX('addtext', 'MEG headshape and sensors',  'Position', [0.5, 0.5, 10,1]);
%exportToPPTX('addpicture', h,  'Position', [1,1,7,5]);

%% Check coregistration (used MNE coregistration, see before)
cfg = [];
mri_orig = ft_determine_coordsys(mri_orig, 'interactive', 'no');
hold on; % add the subsequent objects to the same figure
ft_plot_headshape(shape);
view ([90 0])
title('MEG coregistration', 'FontSize', 13)
%exportToPPTX('addslide');
%exportToPPTX('addpicture', gcf,  'Position',  [1,1,7,5]);

filename = strcat(savemrito, subj, '_mri_orig.mat');
save (filename, 'mri_orig')  % realigned and resampled (if meeded) MRI

%% reslice
mri_orig = ft_convert_units(mri_orig, 'mm');
cfg            = [];
cfg.resolution = 1;
cfg.dim        = [280, 280, 280]; %[256 256 256];
mri_orig_rs    = ft_volumereslice(cfg, mri_orig);
mri_orig = ft_convert_units (mri_orig, 'cm');
mri_orig_rs = ft_convert_units(mri_orig_rs, 'cm');
%mri_orig_realigned = mri_orig_rs;
mri_orig_realigned = mri_orig_rs;
filename = strcat(savemrito, subj, '_mri_orig_realigned.mat');
save (filename, 'mri_orig_realigned')  % realigned and resampled (if meeded) MRI

%% segment
cfg           = [];
cfg.brainthreshold = brainthreshold; % default is 0.5
cfg.output    =   {'brain','skull','scalp'}; %{'brain'};
mri_segmented = ft_volumesegment(cfg, mri_orig_realigned);  
%mri_segmented = ft_determine_coordsys(mri_segmented);
filename = strcat(savemrito, subj, '_mri_segmented.mat');
save (filename, 'mri_segmented')

%% Checking your segmented volumes
%Flipping dimensions after segmenting the volumes (gray, white, and CSF) can easily introduce offsets. 
%Make sure they are correctly aligned with the anatomical scan. See the following procedure.
cfg          = [];
segmentedmri = mri_segmented;
%segmentedmri = ft_volumesegment(cfg, mri); % by default segment white, gray, csf
segmentedmri.transform = mri_orig_realigned.transform;
segmentedmri.anatomy   = mri_orig_realigned.anatomy;
cfg = [];
cfg.method       = 'ortho'; %'slice'
cfg.funparameter = 'brain';
ft_sourceplot(cfg, segmentedmri); 
title (strcat('Segmented volumes, cfg.brainthreshold= ', num2str(brainthreshold) ) );
%exportToPPTX('addslide');
%exportToPPTX('addpicture', gcf,  'Position', [1,1,7,5]);

%%
bss_i                = mri_segmented;
bss_i.seg            = double(mri_segmented.scalp);         % scalp is logical but seg will contain: 0,1,2,3
bss_i.seg(mri_segmented.skull) = 2;                         % skull is represented by index 2 
bss_i.seg(mri_segmented.brain) = 3;                         % brain is represented by index 3
bss_i.seglabel       = {'scalp','skull','brain'}; % label-order corresponds to index from 1 to 3 

cfg                  = [];
cfg.funparameter     = 'seg';
cfg.location         = 'center';
ft_sourceplot(cfg,bss_i);

% change the colormap 
map = [0 0 0; 1 0 0; 0 1 0; 0 0 1];

colormap(map);
title ('Brain,  scull and scalp')
%exportToPPTX('addslide');
%exportToPPTX('addpicture', gcf,  'Position', [1,1,7,5]);

%% headmodel
cfg = [];
cfg.method='singleshell';
individ_hdm_vol = ft_prepare_headmodel(cfg, mri_segmented);
individ_hdm_vol = ft_convert_units(individ_hdm_vol,'cm');

filename = strcat(savemrito, subj,'_individ_hdm_vol.mat');
save (filename, 'individ_hdm_vol')
disp(individ_hdm_vol)

figure;
ft_plot_mesh(individ_hdm_vol.bnd(1),'facecolor','none'); %brain
view ([90 0]);
title ('Individual singleshell head model volume', 'Fontsize', 13)
%exportToPPTX('addslide');
%exportToPPTX('addpicture', gcf,  'Position', [0.5,0.5,7,5]);

%% show brain, sensors and headshape
h=figure
cfg = [];
ft_plot_vol(individ_hdm_vol, 'facecolor' , 'y')
ft_plot_headshape(shape);
ft_plot_sens(grad, 'style', '*b');
title ( 'Brain, sensors and headshape', 'Fontsize', 13);
%exportToPPTX('addslide');
%exportToPPTX('addpicture', gcf,  'Position', [1,1,7,5]);

view ([90 0])
title ( 'Brain, sensors and headshape', 'Fontsize', 13);
%exportToPPTX('addslide');
%exportToPPTX('addpicture', h,  'Position', [1,1,7,5]); %gcf

%% Create template grid based on the standard head model
% note that I save the template in each subj folder (for my convinience),
% but you may do it  once or you may use an availavlke template grid
vol = load(strcat(fieldtripfolder, '/template/headmodel/standard_singleshell'));
cfg = [];
cfg.grid.xgrid  = -25:step:25;
cfg.grid.ygrid  = -25:step:25;
cfg.grid.zgrid  = -25:step:25;
cfg.grid.unit   = 'cm';
cfg.grid.tight  = 'yes';
cfg.inwardshift = -1.5;
cfg.headmodel        = vol.vol;
template_grid  = ft_prepare_sourcemodel(cfg);

filename = strcat(savemrito, subj, '_template_grid.mat');
save (filename, 'template_grid')
% 
h=figure;
ft_plot_mesh(template_grid.pos(template_grid.inside,:));
hold on
ft_plot_vol(vol.vol,  'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
title ('The template grid based on the standard head model', 'Fontsize',13)
%exportToPPTX('addslide');
%exportToPPTX('addpicture', h,  'Position', [0.5,0.5,7,5]);


%% Inverse-warp the subject specific grid to the atlas based template grid
% For this step the individual volume is required.

cfg                = [];
cfg.mri            = mri_orig_realigned;
cfg.grid.warpmni   = 'yes';
cfg.grid.template  =  template_grid; %template.sourcemodel;
cfg.grid.nonlinear = 'yes';%'no';% Be careful here! Check the figure! Sometimes nonlinear transhormation gives strange results
cfg.grid.resolution = 6;
cfg.headmodel = individ_hdm_vol;
%cfg.headshape = shape; 
%cfg.grid.tight = 'yes'
sourcemodel        = ft_prepare_sourcemodel(cfg); % final warped source model

subj_warped_grid = sourcemodel ;
% save warped grid, do it once
filename = strcat(savemrito, subj, '_subj_warped_grid.mat');
save (filename, 'subj_warped_grid')

%% Plot the final source model together with the individual head model and the sensor array

figure; hold on     % plot all objects in one figure

ft_plot_vol(individ_hdm_vol,  'facecolor', 'cortex', 'edgecolor', 'none');  % hdm is individual head model
alpha 0.4           % make the surface transparent
ft_plot_axes(individ_hdm_vol); 
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));% plot only locations inside the volume
%ft_plot_mesh(sourcemodel.pos, 'vertexcolor', 'r');% 
%ft_plot_mesh(template.sourcemodel.pos(template.sourcemodel.inside,:));% plot only locations inside the volume
%ft_plot_mesh(template.sourcemodel.pos, 'vertexcolor', 'g');% plot only locations inside the volume

ft_plot_sens(grad,'style','*r');% plot the sensor array
view ([90 0])
title ('Final source model together with the individual head model and the sensor array', 'FontSize', 13)
%exportToPPTX('addslide');
%exportToPPTX('addpicture', gcf,  'Position', [0.5,0.5,7,5]);


%% Compute the leadfield
% We first create the leadfield using ft_prepare_leadfield using the individual head model from the previous step, the sensor array and the sourcemodel.
% Computes a forward solution for a dipole in a a volume conductor model. The forward solution is expressed as the leadfield matrix (Nchan*3), 
% where each column corresponds with the potential or field distributions on all sensors for one of the x,y,z-orientations of the dipole.
cfg                 = [];
cfg.channel         = grad.label;% ensure that rejected sensors are not present
cfg.grad            = grad;
cfg.headmodel             = individ_hdm_vol;
cfg.lcmv.reducerank = 2; % default for MEG is 2, for EEG is 3
cfg.grid = sourcemodel; % final warped source model
[grid_MNI_lf] = ft_prepare_leadfield(cfg);
filename = strcat(savemrito, subj, '_grid_MNI_lf.mat');
save (filename, 'grid_MNI_lf')

%%
%PPTXname  = strcat(savemrito, subj, '_', date, '_MRI_report');
%exportToPPTX('saveandclose', PPTXname);

%%%
%%%
% mri_convert --in_type mgz --out_type nii --out_orientation RAS /mnt/.../SUBJECT_FOLDER/mri/brainmask.mgz /mnt/.../brainmask.nii.gz[/bash]
% mri_convert --in_type mgz --out_type nii  /Users/mtw/MEG/Karolinska/freesurfersubjects/Case0137/mri/T1.mgz  /Users/mtw/MEG/Karolinska/MRI_nii_subj/0137/0137_T1.nii.gz
% mri_convert --in_type nii  --out_type mgz --out_orientation RAS /Users/mtw/MEG/Karolinska/MRI_nii_subj/0137/x1_0sT1W_3D_32ch_192Hz_mm.nii /Users/mtw/MEG/Karolinska/MRI_nii_subj/0137/0137_T1_192.mgz