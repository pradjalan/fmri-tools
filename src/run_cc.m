function run_cc(input_dir_base,feat_dir_list,ROI_names,mask_thresholds,atlas,output_dir_base)
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


    if ~exist(output_dir_base,'dir')
    [s,c] = system(['mkdir -p ',output_dir_base]);
    end
    
%     if ~isempty(log_path)
%     logfile = fopen(log_path,'w');
%     else
%     logfile = fopen([input_dir_base,'/cc_logs.txt'],'w');
%     end
%     
    if isempty(atlas)
        atlas = ['${FSLDIR}//data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr50-2mm.nii.gz'];
    end
    
    
    
%     [st,num_lines] = system(['cat ' feat_dir_list ' | wc -l ' ]);
%     num_lines = str2num(num_lines);
    
%     fid = fopen(feat_dir_list);
%     fline = fgetl(fid);
    
    flines = strsplit(fileread(feat_dir_list), '\n' );
    num_lines = length(flines);
    
    for cur_line=1:num_lines
%     while ischar(fline)
        fline = flines(cur_line);
        
        scan_loc = strjoin(strcat(input_dir_base,'/',fline));
%         fprintf(logfile,strcat('\nDoing Scan: ',scan_loc));
        disp(scan_loc);
        [~,c] = system(['ls -la ' scan_loc]);
        if isempty(strfind(c,'feat'))==0
            
            feat_loc = strcat(scan_loc,'/.feat/');

            create_transformations(feat_loc);
            
            ROIs = strsplit(ROI_names);
            for rn=1:length(ROIs)
                ROI_name = char(ROIs(rn));
                mask_threshold = mask_thresholds(rn);
                disp(['ROI_NAME: ' ROI_name]);
                disp(['FEAT LOCATION: ' feat_loc]);
%                 fprintf(logfile,'\nExtracting ROIs..');
                
%                 roi_extraction(feat_loc,ROI_name,mask_threshold,atlas,'');
                
                roi_series(feat_loc,ROI_name,mask_threshold,atlas,'');
                
%                 fprintf(logfile,'\nFinding Corelation Maps..');
                %Copy the Directory Structure of Pre-Processed Data
                output_dir = strjoin(strcat(output_dir_base,'/',fline));
                system(['mkdir -p ',output_dir]);
                generate_cc_map(feat_loc,ROI_name,mask_threshold,[1],output_dir);
                
                
                
            end
        end
%         fline = fgetl(fid);
    end
%     fclose(fid);
%     fclose(logfile);

