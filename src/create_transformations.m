function [status]=create_transformations(feat_loc)

% feat_loc : location of feat directory 
% Creates Transformation Matrices for the Subject


% Rotate the highres image so to bring the functional and structural image in the same configration 
 reorient = ['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/fslreorient2std ', feat_loc, '/reg/highres.nii.gz ', feat_loc, '/reg/highres_reoriented.nii.gz "'];
 [status1]=system(reorient);

% Functional Image to Structural Domain. Save the Transformation Matrix for later. 
 Func2High=['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/flirt -in ',feat_loc,'/filtered_func_data.nii.gz',' -ref ',feat_loc,'/reg/highres_reoriented.nii.gz',' -omat ',feat_loc, '/reg/filtered_func2highres_brain.mat"'];
 [status2]=system(Func2High);
 
% Structural Image to Standard Domain. Save Transformation Matrix for later. 
 High2Std=['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/flirt -in ',feat_loc,'/reg/highres_reoriented.nii.gz',' -ref ','${FSLDIR}//data/standard/MNI152_T1_2mm_brain.nii.gz',' -omat ',feat_loc, '/reg/highres_brain2standard.mat"'];
 [status3]=system(High2Std);

% Concatenate Matrices to get Functional to Standard Transformation Matrix.
 ConCat=['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/convert_xfm',' -omat ',feat_loc, '/reg/filtered_func2standard.mat',' -concat ',feat_loc, '/reg/highres_brain2standard.mat ',feat_loc, '/reg/filtered_func2highres_brain.mat"'];
 status4=system(ConCat);

% Invert the matrix to get Standard to Functional Domain Transformation Matrix. 
 Con_xfm=['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/convert_xfm -inverse ',feat_loc, '/reg/filtered_func2standard.mat',' -omat ',feat_loc, '/reg/standard2filtered_func.mat"'];
 status5=system(Con_xfm);

 
% Get Transformation Matrix for Standard to Structural Domain Transformation 
Con_xfm_high=['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/convert_xfm -inverse ',feat_loc, '/reg/highres_brain2standard.mat',' -omat ',feat_loc, '/reg/standard2highres_brain.mat"'];
status9=system(Con_xfm_high);


status=status1+status2+status3+status4+status5+status9;
disp(status);

