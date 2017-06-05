function hypothesis_T_test(input_dir_base,cc_folder_list,ROI_names,output_dir_base)


tom = strsplit('Avg_CC_map_std');
for tn = 1:length(tom)
    type_of_map = char(tom(tn));

    ROIs = strsplit(ROI_names);
    parfor rn=1:length(ROIs)
        ROI_name = char(ROIs(rn));
        disp(['Type of Map: ' type_of_map ';  ROI: ' ROI_name]);
        % type_of_map = Avg_CC_map_std | AvgofMax_CC_map_std | Max_CC_map_std
        fsldir = getenv('FSLDIR');
        file_suffix = [ROI_name];
        out_dir = [output_dir_base,'/',file_suffix,'/'];
        system(['mkdir -p ',out_dir]);
        
        file_name = [out_dir,'/merged_cc_maps_',type_of_map,'.nii.gz'];
%         if ~exist(file_name,'file')
        disp('Creating File list..');
        fid = fopen(cc_folder_list);
        fline = fgetl(fid);
        files = '';
        while ischar(fline)
            fpath = [input_dir_base,'/',fline,'/',file_suffix];
            if exist(fpath,'dir')
            [s,c] = system(['ls -la ',fpath]);
                if exist([fpath '/' type_of_map '.nii.gz'], 'file')
                    files = [files,' ',input_dir_base,'/',fline,'/',file_suffix,'/',type_of_map,'.nii.gz '];
                end
            end
            fline = fgetl(fid);
        end
        % files = fileread(cc_folder_list);
        % files = strrep(files,'\n',' ');
        fclose(fid);
        
        disp('Merging CC_map files..');
        
        merge_cmd = ['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/fslmerge -t ',file_name,' ',files,' "'];
%         disp(merge_cmd);
        
        system(merge_cmd);
        % file_name = [out_dir,'/merged_cc_maps.nii.gz'];
%         end
        
        
        [~,name,extension]=fileparts(file_name);
        out_file_name_T=['T_value_' name extension];
        out_file_name_p=['p_value_' name extension];
        out_file_name_mean=['mean_' name extension];
        out_file_name_stddev=['stddev_' name extension];

        %% Mean and standard deviation of each voxel across all subjects
        disp('Loading Merged Correlation Maps..');
        if ~exist('file_name','file')
            disp(['merged cc map file does not exist... Skipping ROI: ' ROI_name]);
        end
        CC_map_4D=load_untouch_nii(file_name);
        CC_map_4D_img=CC_map_4D.img;

        disp('Creating Brain Mask Variables..');
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
        disp('Calculating Test Statistics of each voxel across subjects..');
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

        
        
        
        disp('Saving Results..');
        % NOTE: Signs are applied to p-values using the t-value matrix
        
        
        pseudo_zero = realmin('double');
        p_CC_map(find(p_CC_map < pseudo_zero)) = pseudo_zero;
        
        nii_T_value_CC_map = std_mask;
        nii_p_value_CC_map = std_mask;
        nii_mean_value_CC_map = std_mask;
        nii_stddev_value_CC_map = std_mask;
        nii_log_p_value_CC_map = std_mask;
        
        nii_T_value_CC_map.hdr=std_mask.hdr; nii_T_value_CC_map.img=T_value_CC_map; save_untouch_nii(nii_T_value_CC_map,[out_dir,'/',out_file_name_T]); %%%% changes to be made in this line
        nii_p_value_CC_map.hdr=std_mask.hdr; nii_p_value_CC_map.img=sign(T_value_CC_map).*abs(p_CC_map); save_untouch_nii(nii_p_value_CC_map,[out_dir,'/',out_file_name_p]);
        disp(['T_value map & p_value map saved in nii format at ',out_dir])

        
        p_CC_map = sign(T_value_CC_map).*(-log10(abs(p_CC_map)));
        nii_log_p_value_CC_map.hdr=std_mask.hdr; nii_log_p_value_CC_map.img=p_CC_map; save_untouch_nii(nii_log_p_value_CC_map,[out_dir,'/log_',out_file_name_p]);
        disp(['log p_value map saved in nii format at ',out_dir])
        
        nii_mean_value_CC_map.hdr=std_mask.hdr; nii_mean_value_CC_map.img=Mean_CC_map; save_untouch_nii(nii_mean_value_CC_map,[out_dir,'/',out_file_name_mean]); %%%% changes to be made in this line
        nii_stddev_value_CC_map.hdr=std_mask.hdr; nii_stddev_value_CC_map.img=Std_dev_CC_map; save_untouch_nii(nii_stddev_value_CC_map,[out_dir,'/',out_file_name_stddev]);
        disp(['Mean map & StdDev map saved in nii format at ',out_dir])


        
        

        % nii_h_CC_map.hdr=std_mask.hdr; nii_h_CC_map.img=h_CC_map; save_nii(nii_h_CC_map,'4D_CC_maps/h_CC_map.nii.gz');
        % nii_p_CC_map.hdr=std_mask.hdr; nii_p_CC_map.img=p_CC_map; save_nii(nii_p_CC_map,'4D_CC_maps/p_CC_map.nii.gz');

    end
end