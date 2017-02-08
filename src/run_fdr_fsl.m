function run_fdr_fsl(input_dir_base,ROI_names,mask_file,output_dir_base)
% Runs the ROI Extraction and Correlation Analysis Code (generate_cc_map) for a given list of feat directories

% input_dir_base : Directory containing the data for all scans

% feat_dir_list : list of locations of feat folders relative to input_dir_base

% ROI_names : names of the ROIs for correlation (matrix object)

% mask_thresholds : atlas intensities of corresponding ROIs

% atlas : Atlas to be used in the analysis (path to nii file). If empty
% string is given, Harvard-Oxford-Cortical (2mm, maxprob50) atlas is used.

% output_dir_base : directory where the correlation results will be saved.
% A new directory structure as per the feat_dir_list will be appended to output_dir_base

% log_path: path to file where the logs for the analysis will be written.
% If empty string is mentioned, logs will be created in input_dir_base. 

% NOTE: A directory is assumed to be present in the feat directory named as
% per the conventions in the previous codes (roi_extraction.m and
% genereate_cc_map.m
    
    if length(output_dir_base)==0
        output_dir_base = input_dir_base;
    end
    
    if ~exist(output_dir_base,'dir')
    [s,c] = system(['mkdir -p ',output_dir_base]);
    end
    
    if length(mask_file)==0
        mask_file = [fsldir,'/data/standard/MNI152_T1_2mm_brain_mask.nii.gz'];    
    end     

    tom = strsplit('Avg_CC_map_std AvgofMax_CC_map_std Max_CC_map_std');
    for tn = 1:length(tom)
        type_of_map = char(tom(tn));
        file_base = ['merged_cc_maps_' type_of_map '.nii.gz'];
        ROIs = strsplit(ROI_names);
            for rn=1:length(ROIs)
                ROI_name = char(ROIs(rn));
                disp(['Type of Map: ' type_of_map ', ROI: ' ROI_name]);
                %Copy the Directory Structure of Pre-Processed Data
                output_dir = strcat(output_dir_base,'/',ROI_name);
                input_dir = strcat(input_dir_base,'/',ROI_name);
                system(['mkdir -p ',output_dir]);
                p_file = [input_dir '/p_value_' file_base];               
                out_file = [output_dir '/fslfdr_' ROI_name '_' type_of_map '.nii.gz'];
                fdr_command = ['fdr -i ' p_file ' -m ' mask_file  ' --conservative -a ' out_file];
                [s,c] = system(fdr_command)
               
                
                
            end
        
        
    end
