

function [Std_CC_Avg_img,Std_CC_Max_img,Std_CC_AvgofMax_img]=generate_cc_map(feat_loc,ROI_name,mask_threshold,type_of_map,output_dir_base)
% feat_loc: Location of feat folder containing functional data and ROI files
% ROI_name: ROI_name, to be appended with output files.
% mask_threshold: Intensity value in the atlas corresponding to ROI
% atlas_threshold: Probability threshold to be used in Harvard-Oxford Atlas
% type_of_map: array containing any of 1,2,3 for specifying type of map to be generated as Average, Max, Avgerage of top 10% respectively
% output_dir_base: Output directory. Output files will be written into a new folder named as per mask_threshold and ROI


disp('Loadind Files for Corelations...');
tic
ROI_dir_name=[feat_loc,'/',ROI_name,'_', num2str(mask_threshold)];
% NOTE: the same suffix is used for reading correlation files in the t-test code.
file_suffix = [ROI_name];
output_dir = [output_dir_base,'/',file_suffix,'/'];
system(['mkdir -p ',output_dir]);

%% Load Input Data

% Load the resting state data
RS=load_untouch_nii([feat_loc,'/filtered_func_data.nii.gz']); %load nifti file of preprocessed functional data
RS_size=size(RS.img);

RS_MASK=load_untouch_nii([feat_loc,'/mask.nii.gz']); %load 3D mask of functional data



% Load ROI data
ROI_MASK=load_untouch_nii([ROI_dir_name,'/ROI_xfmed_mask.nii.gz']); %load transformed mask of ROI
ROI_MASK_img=ROI_MASK.img;
ROI_MASK_size=size(ROI_MASK_img);

ROI_series=load_untouch_nii([ROI_dir_name, '/filtered_func_ROI_masked.nii.gz']); %load masked 4D functional data, time series of the ROI voxels after multiplying with functional data
ROI_size = size(ROI_series.img);



%% Find the coordinates of the ROI and RS mask

disp('Finding Coordinates for ROI and RS..');
% ROI_x=[]; ROI_y=[]; ROI_z=[];
RS_MASK_indices = find(RS_MASK.img);
[RS_x RS_y RS_z] = ind2sub(size(RS_MASK.img), RS_MASK_indices);
old_RS_data = reshape(double(RS.img),RS_size(1)*RS_size(2)*RS_size(3),RS_size(4));
RS_data = zeros(length(RS_MASK_indices),RS_size(4));
RS_data = old_RS_data(RS_MASK_indices,:);


ROI_MASK_indices = find(ROI_MASK.img);
[ROI_x ROI_y ROI_z] = ind2sub(size(ROI_series.img), ROI_MASK_indices);
old_ROI_data = reshape(double(ROI_series.img),ROI_size(1)*ROI_size(2)*ROI_size(3),ROI_size(4));
ROI_data = zeros(length(ROI_MASK_indices),ROI_size(4));
ROI_data = old_ROI_data(ROI_MASK_indices,:);





%normalise
disp('Normalising Data..');
RS_data = zscore(RS_data,0,2);
ROI_data = zscore(ROI_data,0,2);

%Remove Global Signal
disp(' Removing the Global Signal..')
global_signal = nanmean(RS_data,1);
global_component = (global_signal*RS_data.')/(global_signal*global_signal.') ;
RS_data = RS_data - global_component.'*global_signal;
ROI_data = ROI_data - global_component.'*global_signal;



%normalise again
RS_data = zscore(RS_data,0,2);
ROI_data = zscore(ROI_data,0,2);

% % Create Time Series Matrix for ROI
% ROI_indices = find(ROI_MASK_img);
% [ROI_x ROI_y ROI_z] =  ind2sub(ROI_MASK_size, ROI_indices);
% ROI_data = zeros(length(ROI_indices),RS_size(4));
% for i=1:length(ROI_indices)
%     ROI_data(i,:) = RS_data(ROI_indices(i),:);
% end


if isempty(ROI_indices)==0 %i.e. if not empty
    
%% correlation coefficient (CC) of each ROI voxel correspong to each RS voxel
disp('Finding Corelations..');








CORR = (ROI_data*(RS_data.')) / (RS_size(4)-1) ;

%%if u wish to save correlation map in .mat format then uncomment next line
% save([feat_loc,'/',ROI_name,'_',num2str(mask_threshold),'/Correlation_result','CORR','-double'])

%% Avg and Max of Correlation Coef.'s (CC) and of Squareroot of CC's of all ROI voxel's and their mapping
disp('Reducing ROI Corelation Values..');

 %for avg map
if isempty(find(type_of_map==1,1))==0
    disp('Finding Avg ROI Corelation..');
    Avg_CC=nanmean(CORR);%mean(CORR,1,'omitnan');
    Avg_CC_MAP=zeros(size(RS_MASK.img));
    for i=1:length(Avg_CC(1,:))
       Avg_CC_MAP(RS_x(i),RS_y(i),RS_z(i))=Avg_CC(i);
    end
    disp('Wiriting Avg Corelation Files..');
    Avg_CC_img.hdr=RS_MASK.hdr;
    Avg_CC_img.img=Avg_CC_MAP;
    save_nii(Avg_CC_img,[output_dir,'/Avg_CC_map.nii.gz']);
    % registration of map to standard space
    avgmap2std=['sh -c ". ${FSLDIR}//etc/fslconf/fsl.sh;${FSLDIR}/bin/flirt -in ',output_dir,'/Avg_CC_map.nii.gz',' -ref ','${FSLDIR}//data/standard/MNI152_T1_2mm_brain.nii.gz',' -applyxfm -init ',feat_loc, '/reg/filtered_func2standard.mat',' -out ',output_dir,'/Avg_CC_map_std.nii.gz',' -omat ',output_dir,'/Avg_CC_map_2_std.mat"'];
    system(avgmap2std);
    
    func_mask2std=['sh -c ". ${FSLDIR}//etc/fslconf/fsl.sh;${FSLDIR}/bin/flirt -in ',feat_loc,'/mask.nii.gz',' -ref ','${FSLDIR}//data/standard/MNI152_T1_2mm_brain.nii.gz',' -applyxfm -init ',output_dir,'/Avg_CC_map_2_std.mat',' -out ',output_dir,'/MASK_Avg_CC_map_std.nii.gz',' -omat ',output_dir,'/MASK_Avg_CC_map_2_std.mat"'];
    system(func_mask2std);
    %load standard space data of subject to construct 4D map
    Std_CC_Avg=load_untouch_nii([output_dir,'/Avg_CC_map_std.nii.gz']);
    Std_CC_Avg_img=Std_CC_Avg.img;
else
    Std_CC_Avg_img=[];
end

% for max map
if isempty(find(type_of_map==2,1))==0
    disp('Finding Max ROI Corelation..');    
    Max_CC=nanmax(CORR);%max(CORR,[],1,'omitnan');
    Max_CC_neg = nanmax(-CORR);
    Max_CC = sign(Max_CC + Max_CC_neg).*nanmax(Max_CC, abs(Max_CC_neg));
    Max_CC_MAP=zeros(size(RS_MASK.img));
    for i=1:length(Max_CC(1,:))
       Max_CC_MAP(RS_x(i),RS_y(i),RS_z(i))=Max_CC(i);
    end
    
    disp('Writing Max Corelation Files..');
    Max_CC_img.hdr=RS_MASK.hdr;
    Max_CC_img.img=Max_CC_MAP;
    save_nii(Max_CC_img,[output_dir,'/Max_CC_map.nii.gz']);
    
    % registration of map to standard space
    maxmap2std=['sh -c ". ${FSLDIR}//etc/fslconf/fsl.sh;${FSLDIR}/bin/flirt -in ',output_dir,'/Max_CC_map.nii.gz',' -ref ','${FSLDIR}//data/standard/MNI152_T1_2mm_brain.nii.gz',' -applyxfm -init ',feat_loc, '/reg/filtered_func2standard.mat',' -out ',output_dir,'/Max_CC_map_std.nii.gz',' -omat ',output_dir,'/Max_CC_map_2_std.mat"'];
    system(maxmap2std);
    
    func_mask2std=['sh -c ". ${FSLDIR}//etc/fslconf/fsl.sh;${FSLDIR}/bin/flirt -in ',feat_loc,'/mask.nii.gz',' -ref ','${FSLDIR}//data/standard/MNI152_T1_2mm_brain.nii.gz',' -applyxfm -init ',output_dir,'/Max_CC_map_2_std.mat',' -out ',output_dir,'/MASK_Max_CC_map_std.nii.gz',' -omat ',output_dir,'/MASK_Max_CC_map_2_std.mat"'];
    system(func_mask2std);
    %load standard space data of subject to construct 4D map    
    Std_CC_Max=load_untouch_nii([output_dir,'/Max_CC_map_std.nii.gz']);
    Std_CC_Max_img=Std_CC_Max.img;
else
    Std_CC_Max_img=[];
end

% for avgofmax map
if isempty(find(type_of_map==3,1))==0
    disp('Finding AvgOfMax ROI Corelation..');
    Sort_CORR=sort(CORR,1,'descend');
    AvgofMax_CC=mean(Sort_CORR(1:floor((length(CORR(:,1)))/10),:),1); %Average of maximum ~10% correlation coeficient values
    AvgofMax_CC_neg = mean(Sort_CORR(floor((length(CORR(:,1)))*0.9):length(CORR(:,1)),:),1); %Average of maximum ~10% correlation coeficient values
    AvgofMax_CC = sign(AvgofMax_CC + AvgofMax_CC_neg).*nanmax(AvgofMax_CC,abs(AvgofMax_CC_neg));
    AvgofMax_CC_MAP=zeros(size(RS_MASK.img));
    for i=1:length(AvgofMax_CC(1,:))
      AvgofMax_CC_MAP(RS_x(i),RS_y(i),RS_z(i))=AvgofMax_CC(i);
    end
    
    disp('Writing AvgOfMax Corelation Files..');
    AvgofMax_CC_img.hdr=RS_MASK.hdr;
    AvgofMax_CC_img.img=AvgofMax_CC_MAP;
    save_nii(AvgofMax_CC_img,[output_dir,'/AvgofMax_CC_map.nii.gz']);
    
    % registration of map to standard space
    avgmap2std=['sh -c ". ${FSLDIR}//etc/fslconf/fsl.sh;${FSLDIR}/bin/flirt -in ',output_dir,'/AvgofMax_CC_map.nii.gz',' -ref ','${FSLDIR}//data/standard/MNI152_T1_2mm_brain.nii.gz',' -applyxfm -init ',feat_loc, '/reg/filtered_func2standard.mat',' -out ',output_dir,'/AvgofMax_CC_map_std.nii.gz',' -omat ',output_dir,'/AvgofMax_CC_map_2_std.mat"'];
    system(avgmap2std);
    
    func_mask2std=['sh -c ". ${FSLDIR}//etc/fslconf/fsl.sh;${FSLDIR}/bin/flirt -in ',feat_loc,'/mask.nii.gz',' -ref ','${FSLDIR}//data/standard/MNI152_T1_2mm_brain.nii.gz',' -applyxfm -init ',output_dir,'/AvgofMax_CC_map_2_std.mat',' -out ',output_dir,'/MASK_AvgofMax_CC_map_std.nii.gz',' -omat ',output_dir,'/MASK_AvgofMax_CC_map_2_std.mat"'];
    system(func_mask2std);
    %load standard space data of subject to construct 4D map    
    Std_CC_AvgofMax=load_untouch_nii([output_dir,'/AvgofMax_CC_map_std.nii.gz']);
    Std_CC_AvgofMax_img=Std_CC_AvgofMax.img;
else
    Std_CC_AvgofMax_img=[];
end
else
    disp('inappropriate ROI or mask threshold value for ');
    disp('function operation is teminated');
    
end


%% 
toc