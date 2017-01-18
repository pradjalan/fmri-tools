function hypothesis_T_test(file_name,out_dir)
[~,name,extension]=fileparts(file_name);
out_file_name_T=['T_value_' name extension];
out_file_name_p=['p_value_' name extension];
%% Mean and standard deviation of each voxel across all subjects
CC_map_4D=load_untouch_nii(file_name);
CC_map_4D_img=CC_map_4D.img;

std_mask=load_untouch_nii('/usr/local/fsl/data/standard/MNI152_T1_2mm_brain_mask.nii.gz');
std_mask_img=std_mask.img; 
std_mask.hdr.dime.datatype=16; % change header file to use this header for saving T_value in nii format
std_mask.hdr.dime.bitpix=32;
std_mask.hdr.dime.cal_max=0;

mask_indices=find(std_mask_img);
dim_mask_img=[length(std_mask.img(:,1,1)) length(std_mask.img(1,:,1)) length(std_mask.img(1,1,:))]; 
[mask_x,mask_y,mask_z]=ind2sub(dim_mask_img,mask_indices); 
mask_xyz=[mask_x mask_y mask_z]; 

%% Hypothesis test
Mean_CC_map=zeros(dim_mask_img); Std_dev_CC_map=zeros(dim_mask_img); T_value_CC_map=zeros(dim_mask_img);
p_CC_map=zeros(dim_mask_img); %ci_CC_map=zeros(dim_mask_img); stats_CC_map=zeros(dim_mask_img); h_CC_map=zeros(dim_mask_img);


%%
for i=1:length(mask_xyz(:,1))

    Mean_CC_map(mask_x(i,1),mask_y(i,1),mask_z(i,1))=nanmean(CC_map_4D_img(mask_x(i,1),mask_y(i,1),mask_z(i,1),:));
    Std_dev_CC_map(mask_x(i,1),mask_y(i,1),mask_z(i,1))=nanstd(CC_map_4D_img(mask_x(i,1),mask_y(i,1),mask_z(i,1),:));
    T_value_CC_map(mask_x(i,1),mask_y(i,1),mask_z(i,1))=(Mean_CC_map(mask_x(i,1),mask_y(i,1),mask_z(i,1)))/(Std_dev_CC_map(mask_x(i,1),mask_y(i,1),mask_z(i,1))/sqrt(length(CC_map_4D_img(1,1,1,:)))); 

    [~,p_CC_map(mask_x(i,1),mask_y(i,1),mask_z(i,1))] = ttest(CC_map_4D_img(mask_x(i,1),mask_y(i,1),mask_z(i,1),:),0,'Alpha',0.001); %ttest(data,mean,'alpha',alpha_value)
end

% (sample_mean-0)/(std_dev/sqrt(no. of subjects)) [it will give t values of each voxels across all subjects]
%% save file

nii_T_value_CC_map.hdr=std_mask.hdr; nii_T_value_CC_map.img=T_value_CC_map; save_nii(nii_T_value_CC_map,[out_dir,'/',out_file_name_T]); %%%% changes to be made in this line
nii_p_value_CC_map.hdr=std_mask.hdr; nii_p_value_CC_map.img=p_CC_map; save_nii(nii_p_value_CC_map,[out_dir,'/',out_file_name_p]);
disp(['T_value map & p_value map saved in nii format at ',out_dir])
% nii_h_CC_map.hdr=std_mask.hdr; nii_h_CC_map.img=h_CC_map; save_nii(nii_h_CC_map,'4D_CC_maps/h_CC_map.nii.gz');
% nii_p_CC_map.hdr=std_mask.hdr; nii_p_CC_map.img=p_CC_map; save_nii(nii_p_CC_map,'4D_CC_maps/p_CC_map.nii.gz');

