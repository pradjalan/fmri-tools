% addpath('~/Public_writable/MATLAB/NIfTI_20140122/');
addpath('/mnt/project1/rawData/fMRI/incoming/BrainScape_fBIRN/NKI-RS-Lite/COINS/Code/Dixit');
fsldir = '/usr/local/fsl/';
setenv('FSLDIR',fsldir);

fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath); 


filelists = '/mnt/project1/rawData/fMRI/incoming/BrainScape_fBIRN/NKI-RS-Lite/COINS/Filelists/';
logfolder = '/mnt/project1/rawData/fMRI/incoming/BrainScape_fBIRN/NKI-RS-Lite/COINS/CodeResults/logs/';
preprocessed_data_dir = '/mnt/project1/preProcessedData/COINS/withSliceTimingsFile/';
htest_dir = '/mnt/project1/Htests/';
ccmap_dir = '/mnt/project1/CC_maps/COINS/';

atlas_sub = ['${FSLDIR}//data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr50-2mm.nii.gz'];
atlas_cort = ['${FSLDIR}//data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr50-2mm.nii.gz'];
empty_dir = '/mnt/project1/rawData/fMRI/incoming/BrainScape_fBIRN/NKI-RS-Lite/COINS/empty/';
rois = 'mPFC PCC RightAmygdala LeftAmygdala LeftAccumbens RightAccumbens';


% bet_and_flirt_3('/mnt/project1/rawData/fMRI/incoming/BrainScape_fBIRN/NKI-RS-Lite/COINS/empty/.feat8','PCC',30);
% generate_cc_map_3('/mnt/project1/rawData/fMRI/incoming/BrainScape_fBIRN/NKI-RS-Lite/COINS/empty/.feat8','PCC',30,[1,2,3],'/mnt/project1/rawData/fMRI/incoming/BrainScape_fBIRN/NKI-RS-Lite/COINS/empty/CC_Map_check8');
% run_cc('/mnt/project1/preProcessedData/COINS/withSliceTimingsFile/','/mnt/project1/rawData/fMRI/incoming/BrainScape_fBIRN/NKI-RS-Lite/COINS/Filelists/healthy_subjects_adults_scan_directories_1400.txt','PCC',30,'/mnt/project1/CC_maps/COINS/DMN_all/preProcessedwithSliceTimingsFile/','/mnt/project1/rawData/fMRI/incoming/BrainScape_fBIRN/NKI-RS-Lite/COINS/CodeResults/logs/cc_map/DMN_all_healthy/2/dmn_all_healthy_1400.txt')