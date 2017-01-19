function fdr_correction(t_fileLoc,p_fileLoc,q_value,out_dir)
% fdr_correction(p_fileLoc,q_value,out_dir)
%%FDR correction
fsldir = getenv('FSLDIR');
[~,name,extension]=fileparts(p_fileLoc);
out_file_name=['fdr_corr_' name extension]; % specify output file name

p_map_nii=load_untouch_nii(p_fileLoc);
p_map_img=p_map_nii.img; %seperating img file from nii file

std_brain_mask=load_untouch_nii([fsldir '/data/standard/MNI152_T1_2mm_brain_mask.nii.gz']);
std_mask_img=double(std_brain_mask.img); %seperating img file from standard nii file

mask_indices=find(std_mask_img); %to find indices of brain mask

p_values=p_map_img(mask_indices); %extracting p_values of brain only (i.e. excluding outside brain region)
% bcz p value of outside brain location is also zero
p_val_mask=p_values; p_val_mask(:)=1; %p value mask is required when to align sorted matrix back to its original position
[sorted_p_values,sorted_p_indices]=sort(p_values,'ascend');

mafd = mafdr(p_values);
LENGTH=length(sorted_p_values);

for i=1:length(sorted_p_values)
    if sorted_p_values(i)<(i/LENGTH)*q_value
        continue
    end
    r=i
    sorted_p_values(i)
    LENGTH
    p_val_mapped_on_brain=p_val_mask; p_val_mapped_on_brain(sorted_p_indices(1:r))=sorted_p_values(1:r); 
    break
end
        %%
std_mask_img(mask_indices)=p_val_mapped_on_brain;
fdr_corr_nii.img=std_mask_img; fdr_corr_nii.hdr=p_map_nii.hdr;
save_nii(fdr_corr_nii,[out_dir,'/',out_file_name])

std_mask_img(mask_indices)=mafd;
mat_fdr_corr_nii.img=std_mask_img; mat_fdr_corr_nii.hdr=p_map_nii.hdr;
save_nii(mat_fdr_corr_nii,[out_dir,'/mat_',out_file_name])




CreateMask=['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/fslmaths ',out_dir,'/',out_file_name,' -thr 0.5 -binv ',out_dir,'/binv_',out_file_name,' "'];
scm = system(CreateMask)

ApplyMask = ['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/fslmaths ',t_fileLoc,' -mul ',out_dir,'/binv_',out_file_name,'  ',out_dir,'/fdr_corrected_',out_file_name,'  "'];
sam = system(ApplyMask)

disp(['FDR corrected file with name ',out_file_name,' is saved at ', out_dir])
