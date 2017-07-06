function compute_roi_intersection(input_dir,rois,rois_dir,output_dir,mask_file)
%config file contains coordinates, radii and names


fsldir = getenv('FSLDIR');

if length(mask_file)==0
    mask_file = [fsldir,'/data/standard/MNI152_T1_2mm_brain_mask.nii.gz'];
end


% rois = strrep(strsplit(strtrim(ls(rois_dir))),'.nii.gz','');
rois = strsplit(rois);

% lsr = 'AccumbensR_Addiction_3.nii.gz AccumensL_Addiction_3.nii.gz AmygdalaL_Addiction_4.nii.gz AmygdalaR_Addiction_4.nii.gz HippocampusL_Addiction_4.nii.gz HippocampusR_Addiction_4.nii.gz IPLL_DMN_4.nii.gz IPLR_DMN_4.nii.gz MPFC_DMN_3.nii.gz PCC_HO_4.nii.gz';
% rois = strrep(strsplit(strtrim(lsr)),'.nii.gz','');

parameters = strsplit('TotalVoxels,NonZeroVoxels,Volume,Mean,Percentile100(Max),Percentile90,Percentile75,Percentile50(Median),Percentile25,Percentile10,Percentile0(Min),MeanPositive,MeanNegative',',');

out_matrix = zeros(length(rois),length(rois), length(parameters));

for r=1:length(rois)
    cur_roi = strrep(char(rois(r)),'.nii.gz','');
    disp(cur_roi);
%     input_file = [input_dir '/' cur_roi '_logq_value0002.nii.gz'];
    input_file = [input_dir '/' cur_roi '/q_values_' cur_roi '_Avg_CC_map_std.nii.gz'];
% fid = fopen(out_file,'w');
% fprintf(fid,'ROI_name,TotalVoxels,NonZeroVoxels,Volume,Mean,Percentile100(Max),Percentile90,Percentile75,Percentile50(Median),Percentile25,Percentile10,Percentile0(Min),MeanPositive,MeanNegative\n');


    for rn=1:length(rois)
     roi = strrep(char(rois(rn)),'.nii.gz','');
     roi_mask_path = [rois_dir '/' roi '.nii.gz'];

     stat_command = ['fslstats ' roi_mask_path ' -V' ];
     [s,c] = system(['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/' stat_command '  "']);
     if s~=0
        disp('Error in fslstats..');
     end

     roi_voxels = strsplit(c);
     roi_voxels = c(1);

     mul_command = ['fslmaths ' input_file ' -mul ' roi_mask_path ' ' rois_dir '/temp_blob_intersection.nii.gz' ];
     mul_neg_command = ['fslmaths ' input_file ' -uthr 0 -mul ' roi_mask_path ' ' rois_dir '/temp_blob_intersection_negative.nii.gz' ];
     mul_pos_command = ['fslmaths ' input_file ' -thr 0 -mul ' roi_mask_path ' ' rois_dir '/temp_blob_intersection_positive.nii.gz' ];

     
     [s,~] =  system(['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/' mul_command '  "']);
        if s~=0
            disp('Problem Creating ROI series:');
            disp(mul_command);
        end

    [s,~] =  system(['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/' mul_neg_command '  "']);
        if s~=0
            disp('Problem Creating Negative ROI series:');
        end

    [s,~] =  system(['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/' mul_pos_command '  "']);
        if s~=0
            disp('Problem Creating Positive ROI series:');
        end

        disp(['ROI: ' roi ';  ']);




     stat_command = ['fslstats ' rois_dir '/temp_blob_intersection.nii.gz' ' -V -M -P 100 -P 90 -P 75 -P 50 -P 25 -P 10 -P 0' ];
     [s,c] = system(['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/' stat_command '  "']);
     if s~=0
        disp('Error Running fslstats..');
     end
     file_line = strjoin(strsplit(c),',');
     mat_line = str2num([roi_voxels ',' file_line(1:end-1)]);
     
     

     stat_command_pos = ['fslstats ' rois_dir '/temp_blob_intersection_positive.nii.gz' ' -M' ];
     [s,c] = system(['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/' stat_command_pos '  "']);
     if s~=0
        disp('Error Running fslstats..');
     end
     mat_line = [mat_line  str2num(strtrim(c))];



     stat_command_neg = ['fslstats ' rois_dir '/temp_blob_intersection_negative.nii.gz' ' -M' ];
     [s,c] = system(['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/' stat_command_neg '  "']);
     if s~=0
        disp('Error Running fslstats..');
     end
     mat_line = [mat_line  str2num(strtrim(c))];
     

     
     
     out_matrix(r,rn,:) = mat_line;
%      fprintf(fid,[file_line '\n']);
     delete([rois_dir '/temp_blob_intersection*']);
     
    end
end
% fclose(fid);     

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

save([output_dir '/data_matrix.mat'],'out_matrix'); 


%%%%%%%%%%%%%% USE imagesc(matrix) to plot 2D matrices...

for p = 1:length(parameters)
    figure('visible','off');
    M = out_matrix(:,:,p);
    imagesc(M.', [-10 10]);
    set(gca,'Ytick',1:length(rois));
    set(gca,'YtickLabel',strrep(rois,'_','\_'));
    title(parameters(p));
    colorbar;
    saveas(gca,[output_dir '/' char(parameters(p)) '.png']);
end

end