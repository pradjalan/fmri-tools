function two_sided_t_test(group1_rootdir,group2_rootdir,ROI_names,output_dir_base)

types_of_map = strsplit('Avg_CC_map_std AvgofMax_CC_map_std Max_CC_map_std');
ROIs = strsplit(ROI_names);
for tn=1:length(types_of_map)
    type_of_map = char(types_of_map(tn));
    for rn=1:length(ROIs)
        ROI_name = char(ROIs(rn));
        disp(ROI_name);
        disp(type_of_map);
        fsldir = getenv('FSLDIR');
        file_suffix = [ROI_name];
        group1_dir = [group1_rootdir,'/',file_suffix,'/'];
        group2_dir = [group2_rootdir,'/',file_suffix,'/'];
        out_dir = [output_dir_base '/' file_suffix '/'];
        disp(out_dir);
        system(['mkdir -p ' out_dir]);
        
        
        
        group1_mean_file_name = [group1_dir,'/mean_merged_cc_maps_',type_of_map,'.nii.gz'];
        group2_mean_file_name = [group2_dir,'/mean_merged_cc_maps_',type_of_map,'.nii.gz'];
        group1_stddev_file_name = [group1_dir,'/stddev_merged_cc_maps_',type_of_map,'.nii.gz'];
        group2_stddev_file_name = [group2_dir,'/stddev_merged_cc_maps_',type_of_map,'.nii.gz'];
        
        cc_group1_file = [group1_dir,'/merged_cc_maps_',type_of_map,'.nii.gz'];
        cc_group2_file = [group2_dir,'/merged_cc_maps_',type_of_map,'.nii.gz'];
        disp(cc_group1_file);
        
        
        [~,name,extension]=fileparts(cc_group1_file);
        out_file_name_T=['T_value_' type_of_map '.nii.gz']
        out_file_name_p=['p_value_' type_of_map '.nii.gz'];
%         out_file_name_mean=['mean_' name extension];
%         out_file_name_stddev=['stddev_' name extension];

        %% Mean and standard deviation of each voxel across all subjects
        hm_map=load_untouch_nii(group1_mean_file_name);
        hs_map=load_untouch_nii(group1_stddev_file_name);
        am_map=load_untouch_nii(group2_mean_file_name);
        as_map=load_untouch_nii(group2_stddev_file_name);
        
        CC_map_4D_group1=load_untouch_nii(cc_group1_file);
        CC_map_4D_group1_img=CC_map_4D_group1.img;
        
        CC_map_4D_group2=load_untouch_nii(cc_group2_file);
        CC_map_4D_group2_img=CC_map_4D_group2.img;
        
        
        std_mask=load_untouch_nii([fsldir,'/data/standard/MNI152_T1_2mm_brain_mask.nii.gz']);
        std_mask_img=std_mask.img; 
        std_mask.hdr.dime.datatype=16; % change header file to use this header for saving T_value in nii format
        std_mask.hdr.dime.bitpix=32;
        std_mask.hdr.dime.cal_max=0;

        mask_indices=find(std_mask_img);
        dim_mask_img=[length(std_mask.img(:,1,1)) length(std_mask.img(1,:,1)) length(std_mask.img(1,1,:))]; 
        [mask_x,mask_y,mask_z]=ind2sub(dim_mask_img,mask_indices); 
        mask_xyz=[mask_x mask_y mask_z]; 

        %% Hypothesis test
        Mean_CC_map=zeros(dim_mask_img); 
        Std_dev_CC_map=zeros(dim_mask_img); 
        T_value_CC_map=zeros(dim_mask_img);
        p_CC_map=zeros(dim_mask_img); %ci_CC_map=zeros(dim_mask_img); stats_CC_map=zeros(dim_mask_img); h_CC_map=zeros(dim_mask_img);

        
        %%
        for i=1:length(mask_xyz(:,1))

%             Mean_CC_map(mask_x(i,1),mask_y(i,1),mask_z(i,1))=nanmean(CC_map_4D_img(mask_x(i,1),mask_y(i,1),mask_z(i,1),:));
%             Std_dev_CC_map(mask_x(i,1),mask_y(i,1),mask_z(i,1))=nanstd(CC_map_4D_img(mask_x(i,1),mask_y(i,1),mask_z(i,1),:));
            
            nume = (hm_map.img(mask_x(i,1),mask_y(i,1),mask_z(i,1)) - am_map.img(mask_x(i,1),mask_y(i,1),mask_z(i,1)) );
            n1 = length(CC_map_4D_group1_img(1,1,1,:));
            n2 = length(CC_map_4D_group2_img(1,1,1,:));
            deno = sqrt( ((hs_map.img(mask_x(i,1),mask_y(i,1),mask_z(i,1))).^2)/(n1)  +  ((hs_map.img(mask_x(i,1),mask_y(i,1),mask_z(i,1))).^2)/(n2) ); 
            T_value_CC_map(mask_x(i,1),mask_y(i,1),mask_z(i,1)) = nume/deno;
            
            [~,p_CC_map(mask_x(i,1),mask_y(i,1),mask_z(i,1))] = ttest2(CC_map_4D_group1_img(mask_x(i,1),mask_y(i,1),mask_z(i,1),:), CC_map_4D_group2_img(mask_x(i,1),mask_y(i,1),mask_z(i,1),:),'Alpha',0.001); %ttest(data,mean,'alpha',alpha_value)
            
%             [~,p_CC_map(mask_x(i,1),mask_y(i,1),mask_z(i,1))] = ttest(CC_map_4D_img(mask_x(i,1),mask_y(i,1),mask_z(i,1),:),0,'Alpha',0.001); %ttest(data,mean,'alpha',alpha_value)
%             length(CC_map_4D_group1_img(mask_x(i,1),mask_y(i,1),mask_z(i,1),:))
%             length(CC_map_4D_group2_img(mask_x(i,1),mask_y(i,1),mask_z(i,1),:))
%             [h,p_CC_map(mask_x(i,1),mask_y(i,1),mask_z(i,1)),ci,stats] = ttest(CC_map_4D_group1_img(mask_x(i,1),mask_y(i,1),mask_z(i,1),:), CC_map_4D_group2_img(mask_x(i,1),mask_y(i,1),mask_z(i,1),:) , 'Alpha',0.001);
%             T_value_CC_map(mask_x(i,1),mask_y(i,1),mask_z(i,1)) = stats.tstat;
        end

        % (sample_mean-0)/(std_dev/sqrt(no. of subjects)) [it will give t values of each voxels across all subjects]
        %% save file
        disp([out_dir,'/',out_file_name_T]);
        nii_T_value_CC_map.hdr=std_mask.hdr; nii_T_value_CC_map.img=T_value_CC_map; save_nii(nii_T_value_CC_map,[out_dir,'/',out_file_name_T]); %%%% changes to be made in this line
        nii_p_value_CC_map.hdr=std_mask.hdr; nii_p_value_CC_map.img=p_CC_map; save_nii(nii_p_value_CC_map,[out_dir,'/',out_file_name_p]);
        disp(['T_value and p_value map saved in nii format at: ',out_dir])

        pseudo_zero = realmin('double');
        p_CC_map(find(p_CC_map < pseudo_zero)) = pseudo_zero;
        p_CC_map = sign(T_value_CC_map).*(-log10(abs(p_CC_map)));
        nii_log_p_value_CC_map.hdr=std_mask.hdr; nii_log_p_value_CC_map.img=p_CC_map; save_nii(nii_log_p_value_CC_map,[out_dir,'/log_',out_file_name_p]);
        disp(['log p_value map saved in nii format at ',out_dir])
        
        
%         nii_mean_value_CC_map.hdr=std_mask.hdr; nii_mean_value_CC_map.img=Mean_CC_map; save_nii(nii_mean_value_CC_map,[out_dir,'/',out_file_name_mean]); %%%% changes to be made in this line
%         nii_stddev_value_CC_map.hdr=std_mask.hdr; nii_stddev_value_CC_map.img=Std_dev_CC_map; save_nii(nii_stddev_value_CC_map,[out_dir,'/',out_file_name_stddev]);
%         disp(['T_value map & p_value map saved in nii format at ',out_dir])



        % nii_h_CC_map.hdr=std_mask.hdr; nii_h_CC_map.img=h_CC_map; save_nii(nii_h_CC_map,'4D_CC_maps/h_CC_map.nii.gz');
        % nii_p_CC_map.hdr=std_mask.hdr; nii_p_CC_map.img=p_CC_map; save_nii(nii_p_CC_map,'4D_CC_maps/p_CC_map.nii.gz');

    end
end