function [status]=roi_series_manual_reg(feat_loc,ROI_name,mask_threshold,atlas,output_dir)

% feat_loc : location of feat directory 
%
% ROI_name : name of ROI (this name will be appended in the output directories for further analysis
%
% mask_threshold : Intensity value assigned to the ROI in the atlas being used. 
%
% atlas : path to the atlas file (in nii format). If empty string is
% mentioned, then Harvard-Oxford-Cortical (2mm) atlas is used.
%
% output_dir: All the output matrices and maps are written to this
% directory. If empty string is mentioned then a directory named using
% ROI_name and mask_threshold is created in feat_loc directory.


if isempty(atlas)
    atlas = ['${FSLDIR}//data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr50-2mm.nii.gz'];
end

if isempty(output_dir)
    ROI_dir_name=[feat_loc,'/rois/',ROI_name,'_',num2str(mask_threshold)];
    mkdir(ROI_dir_name);
else
    ROI_dir_name = output_dir;
end



% Extract the ROI Mask from the given atlas using intensity index (mask_threshold)
 Xfmation=['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/fslmaths ', atlas, ' -thr ',num2str(mask_threshold),' -uthr ', num2str(mask_threshold),' -bin ',ROI_dir_name,'/ROI_mask.nii.gz"'];
 status9=system(Xfmation);

% Bring the ROI Mask to Functional Domain
 FLIRT=['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/flirt -in ',ROI_dir_name,'/ROI_mask.nii.gz',' -ref ',feat_loc,'/filtered_func_data.nii.gz',' -applyxfm -init ',feat_loc, '/reg/standard2filtered_func.mat',' -out ',ROI_dir_name,'/ROI_xfmed.nii.gz"'];
 status6=system(FLIRT);
 
% Remove the effect of smoothening to exclude extra voxels in binary mask
Xfmation=['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/fslmaths ',ROI_dir_name,'/ROI_xfmed.nii.gz -thr ',num2str(0.5),' -bin ',ROI_dir_name,'/ROI_xfmed_mask.nii.gz"'];
status7=system(Xfmation);

% Apply ROI Mask to image in the functional domain
xfm_Mask=['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/fslmaths ',feat_loc,'/filtered_func_data.nii.gz -mul ',ROI_dir_name,'/ROI_xfmed_mask.nii.gz ',ROI_dir_name,'/filtered_func_ROI_masked.nii.gz"'];
status8=system(xfm_Mask);

% Save ROI 
FLIRT_high=['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/flirt -in ',ROI_dir_name,'/ROI_mask.nii.gz',' -ref ',feat_loc,'/reg/highres_reoriented.nii.gz',' -applyxfm -init ',feat_loc,'/reg/standard2highres_brain.mat',' -out ',ROI_dir_name,'/ROI_mask_highres.nii.gz"'];
status10=system(FLIRT_high);

status=status6+status7+status8+status9+status10;
disp(status);

