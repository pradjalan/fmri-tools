function run_cc_parallel(input_dir_base,feat_dir_list,ROIs_dir, ROI_names,output_dir_base)
        ROIs = strsplit(ROI_names);
        for rn=1:length(ROIs)
            ROI_name = char(ROIs(rn));
            roi_atlas = [ROIs_dir '/' ROI_name '.nii.gz' ];
            run_cc(input_dir_base,feat_dir_list,ROI_name,[1],roi_atlas,output_dir_base);
        end
end