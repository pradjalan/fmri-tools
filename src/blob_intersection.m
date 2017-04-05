function blob_intersection(input_file,rois_dir,rois,out_file,mask_file)
%config file contains coordinates, radii and names


fsldir = getenv('FSLDIR');

if length(mask_file)==0
    mask_file = [fsldir,'/data/standard/MNI152_T1_2mm_brain_mask.nii.gz'];
end
% 
% if ~exist(out_dir,'dir')
%     system(['mkdir -p ' out_dir]);
% end

data_file = load_untouch_nii(input_file);
% image_matrix = data_file.img;
% 
% standard_header=load_untouch_header_only(mask_file);
% standard_file = load_untouch_nii(mask_file);
% standard_image = standard_file.img;

rois = strsplit(rois);

fid = fopen(out_file,'w');
fprintf(fid,'ROI_name,NonZeroVoxels,Volume,Mean,Percentile100(Max),Percentile90,Percentile75,Percentile50(Median),Percentile25,Percentile10,Percentile0(Min),MeanPositive,MeanNegative\n');



for rn=1:length(rois)
 roi = char(rois(rn));
 roi_mask_path = [rois_dir '/' roi '_3.nii.gz'];
 mul_command = ['fslmaths ' input_file ' -mul ' roi_mask_path ' ' rois_dir '/temp_blob_intersection.nii.gz' ];
 mul_neg_command = ['fslmaths ' input_file ' -uthr 0 -mul ' roi_mask_path ' ' rois_dir '/temp_blob_intersection_negative.nii.gz' ];
 mul_pos_command = ['fslmaths ' input_file ' -thr 0 -mul ' roi_mask_path ' ' rois_dir '/temp_blob_intersection_positive.nii.gz' ];

 [s,~] =  system(['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/' mul_command '  "']);
    if s~=0
        disp('Problem Creating ROI series:');
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
 
 stat_command_pos = ['fslstats ' rois_dir '/temp_blob_intersection_positive.nii.gz' ' -M' ];
 [s,c] = system(['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/' stat_command_pos '  "']);
 if s~=0
    disp('Error Running fslstats..');
 end
 file_line = [file_line ',' c];



 stat_command_neg = ['fslstats ' rois_dir '/temp_blob_intersection_negative.nii.gz' ' -M' ];
 [s,c] = system(['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/' stat_command_pos '  "']);
 if s~=0
    disp('Error Running fslstats..');
 end
 file_line = [file_line ',' c];
 
 fprintf(fid,[file_line '\n']);
 
 
end

fclose(fid);     


end