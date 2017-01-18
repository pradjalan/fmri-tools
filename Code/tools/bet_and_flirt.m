function [status]=bet_and_flirt(feat_loc,ROI_loc,ROI_name,mask_threshold,bet_threshold)
% bet_and_flirt(feat_loc,ROI_loc,ROI_name,mask_threshold,bet_threshold)
% feat_loc - location of feat directory (including feat directory name 
% e.g.'~/Documents/cerebrum_data/NKI-RS-Lite/Release_5/group_4/0114326/session_1/RfMRI_std_2500.feat')
% ROI_loc - location of ROI mask file 
% ROI_name - name of ROI
% mask_thresh - threshold for binary mask (0-1) or probablistic mask (1-100)

ROI_dir_name=[feat_loc,'/',ROI_name,'_',num2str(mask_threshold)];
mkdir(ROI_dir_name);
atlas = '${FSLDIR}/data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz';
% brain extraction with given threshold value
% if bet_threshold~=0
%     BET=['${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/bet ',feat_loc,'/reg/highres.nii.gz ',feat_loc,'/reg/highres_brain.nii.gz',' -f ',num2str(bet_threshold),'"'];
%     [status1]= system(BET);
% end
status1 = 0
disp(system(['echo $FSLDIR']))

% 

Func2High=['sh -c ". ${FSLDIR}etc/fslconf/fsl.sh;${FSLDIR}bin/flirt -in ',feat_loc,'/filtered_func_data.nii.gz',' -ref ',feat_loc,'/reg/highres.nii.gz',' -omat ',ROI_dir_name,'/filtered_func2highres_brain.mat"'];
[status2]=system(Func2High)

High2Std=['sh -c ". ${FSLDIR}etc/fslconf/fsl.sh;${FSLDIR}bin/flirt -in ',feat_loc,'/reg/highres.nii.gz',' -ref ','${FSLDIR}/data/standard/MNI152_T1_2mm_brain.nii.gz',' -omat ',ROI_dir_name,'/highres_brain2standard.mat"'];
[status3]=system(High2Std)
% system('sh -c ". ${FSLDIR}etc/fslconf/fsl.sh;${FSLDIR}bin/flirt -in ../reg/highres_brain.nii.gz -ref ../reg/standard.nii.gz -omat highres_brain2standard.mat"')

ConCat=['sh -c ". ${FSLDIR}etc/fslconf/fsl.sh;${FSLDIR}bin/convert_xfm',' -omat ',ROI_dir_name,'/filtered_func2standard.mat',' -concat ',ROI_dir_name,'/highres_brain2standard.mat ',ROI_dir_name,'/filtered_func2highres_brain.mat"'];
status4=system(ConCat)

Con_xfm=['sh -c ". ${FSLDIR}etc/fslconf/fsl.sh;${FSLDIR}bin/convert_xfm -inverse ',ROI_dir_name,'/filtered_func2standard.mat',' -omat ',ROI_dir_name,'/standard2filtered_func.mat"'];
status5=system(Con_xfm)

FLIRT=['sh -c ". ${FSLDIR}etc/fslconf/fsl.sh;${FSLDIR}bin/flirt -in ',ROI_loc,'/',ROI_name,' -ref ',feat_loc,'/filtered_func_data.nii.gz',' -applyxfm -init ',ROI_dir_name,'/standard2filtered_func.mat',' -out ',ROI_dir_name,'/ROI_xfmed.nii.gz"'];
status6=system(FLIRT)

Xfmation=['sh -c ". ${FSLDIR}etc/fslconf/fsl.sh;${FSLDIR}bin/fslmaths ',ROI_dir_name,'/ROI_xfmed.nii.gz -thr ',num2str(mask_threshold),' -bin ',ROI_dir_name,'/ROI_xfmed_mask.nii.gz"'];
status7=system(Xfmation)

xfm_Mask=['sh -c ". ${FSLDIR}etc/fslconf/fsl.sh;${FSLDIR}bin/fslmaths ',feat_loc,'/filtered_func_data.nii.gz -mul ',ROI_dir_name,'/ROI_xfmed_mask.nii.gz ',ROI_dir_name,'/filtered_func_ROI_masked.nii.gz"'];
status8=system(xfm_Mask)

%register ROI on highres.nii.gz image
Con_xfm_high=['sh -c ". ${FSLDIR}etc/fslconf/fsl.sh;${FSLDIR}bin/convert_xfm -inverse ',ROI_dir_name,'/highres_brain2standard.mat',' -omat ',ROI_dir_name,'/standard2highres_brain.mat"'];
status9=system(Con_xfm_high)

FLIRT_high=['sh -c ". ${FSLDIR}etc/fslconf/fsl.sh;${FSLDIR}bin/flirt -in ',ROI_loc,'/',ROI_name,' -ref ',feat_loc,'/reg/highres.nii.gz',' -applyxfm -init ',ROI_dir_name,'/standard2highres_brain.mat',' -out ',ROI_dir_name,'/ROI_highres.nii.gz"'];
status10=system(FLIRT_high)
status=status1+status2+status3+status4+status5+status6+status7+status8+status9+status10;

