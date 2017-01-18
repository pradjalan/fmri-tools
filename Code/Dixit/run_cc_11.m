function run_cc_11(input_dir_base,feat_dir_list,ROI_name,mask_threshold,atlas_threshold,output_dir_base,log_path)
% Designed for running on 1000FARMs data
    [s,c]=system(['mkdir -p ',output_dir_base])
%     ROIs = strsplit(ROI_names);
%     for rn = 1:length(ROIs)
%         ROI_name = char(ROIs(rn));
        logfile = fopen(log_path,'w');
        fid = fopen(feat_dir_list);
        fline = fgetl(fid);
        num_scans = 0;
        while ischar(fline)
            scan_loc = strcat(input_dir_base,'/',fline);
            fprintf(logfile,strcat('\nDoing Scan: ',scan_loc));
            [s,c] = system(['ls -la ',scan_loc]);
            if isempty(strfind(c,'feat'))==0
                feat_loc = strcat(scan_loc,'/rest.feat/');
                fprintf(logfile,'\nExtracting ROIs..');
                bet_and_flirt_3_1(feat_loc,ROI_name,mask_threshold,atlas_threshold);
                fprintf(logfile,'\nFinding Corelation Maps..');
                output_dir = strcat(output_dir_base,'/',fline);
                system(['mkdir -p ',output_dir]);
                generate_cc_map_3(feat_loc,ROI_name,mask_threshold,atlas_threshold,[1,2,3],output_dir);
                num_scans = num_scans + 1;
                fprintf(logfile,strcat('\nScans Completed:',num2str(num_scans)) );
            end
            fline = fgetl(fid);
        end
        fclose(fid);
        fclose(logfile);

%     end