function run_cc_parallel(input_dir_base,feat_dir_list,ROI_names,output_dir_base,log_path)
        ROIs = strsplit(ROI_names);
        parfor rn=1:length(ROIs)
            ROI_name = char(ROIs(rn));
            roi_atlas = ['/mnt/project1/home1/cs5120287/data/pc_rois/' ROI_name '.nii.gz' ];
            run_cc(input_dir_base,feat_dir_list,ROI_name,[1],roi_atlas,output_dir_base,log_path);
        end

end