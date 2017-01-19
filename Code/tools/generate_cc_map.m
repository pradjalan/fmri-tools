
function [Std_CC_Avg_img,Std_CC_Max_img,Std_CC_AvgofMax_img]=generate_cc_map(feat_loc,ROI_name,mask_threshold,type_of_map,output_dir)
% Use it after bet_and_flirt function
% All inputs are required & input must be a string
% Corr_map(rs,rs_masked,roi_MASK,rs_MASK)
% rs='../filtered_func_data.nii.gz'
% rs_masked='filtered_func_ROI_masked.nii.gz'
% roi_MASK='ROI_xfmed_MASK.nii.gz'
% rs_MASK='../mask.nii.gz'
% DIR_path=<path of directory> 
disp('generating correlation maps...')
tic
ROI_dir_name=[feat_loc,'/',ROI_name,'_',num2str(mask_threshold)];
% Loc=strcat('Users/dixit/Documents/cerebrum_data/',num2str(
%% Input Data
RS=load_untouch_nii([feat_loc,'/filtered_func_data.nii.gz']); %load nifti file of preprocessed functional data

RS_masked=load_untouch_nii([feat_loc,'/',ROI_name,'_',num2str(mask_threshold),'/filtered_func_ROI_masked.nii.gz']); %load masked 4D functional data

ROI_MASK=load_untouch_nii([feat_loc,'/',ROI_name,'_',num2str(mask_threshold),'/ROI_xfmed_mask.nii.gz']); %load transformed mask of ROI
ROI_MASK_img=ROI_MASK.img;

RS_MASK=load_untouch_nii([feat_loc,'/mask.nii.gz']); %load 3D mask of functional data


%% Find the coordinates of the ROI mask
ROI_x=[]; ROI_y=[]; ROI_z=[];
for i=1:length(ROI_MASK.img(1,1,:))
    clear z1
    [x1,y1]=find(ROI_MASK_img(:,:,i));
    if x1>0
    z1(1:length(x1),1)=i;
    ROI_x=[ROI_x;x1]; ROI_y=[ROI_y;y1]; ROI_z=[ROI_z;z1];
    end
end
clear x1 y1 z1
ROI_MASK_xyz=[ROI_x ROI_y ROI_z];


if isempty(ROI_MASK_xyz)==0 %i.e. if not empty

%%Time series (TS) of voxels of seed region (ROI)
RS_ROI_vox_TS=zeros(1,length(RS.img(1,1,1,:)));
for j=1:length(ROI_MASK_xyz(:,1))
    for  i=1:length(RS.img(1,1,1,:))
    
    RS_ROI_vox_TS(j,i)=RS_masked.img(ROI_x(j),ROI_y(j),ROI_z(j),i); % to find time series of each voxel of seed region
    end
end

%%time series (TS) of each voxel of brain region

% ROI_MASK_indices=find(ROI_MASK.img);
RS_MASK_indices=find(RS_MASK.img);

RS_vox_TS=zeros(length(RS_MASK_indices),length(RS.img(1,1,1,:)));

% conversion of indices value to [x y z] coordinate value for RS_MASK
S=[length(RS_MASK.img(:,1,1)) length(RS_MASK.img(1,:,1)) length(RS_MASK.img(1,1,:))]; %dimensions of RS_MASK.img
[RS_x,RS_y,RS_z]=ind2sub(S,RS_MASK_indices); %conversion of indices to coordinate value
RS_MASK_xyz=[RS_x RS_y RS_z]; %concatenating 


for j=1:length(RS_MASK_xyz(:,1))
    for i=1:length(RS.img(1,1,1,:))
        
        RS_vox_TS(j,i)=RS.img(RS_x(j),RS_y(j),RS_z(j),i);
    end
end

    
%% correlation coefficient (CC) of each ROI voxel correspong to each RS voxel
StdDev_ROI_TS=std(RS_ROI_vox_TS,0,2);
StdDev_RS_TS=std(RS_vox_TS,0,2);
RS_ROI_vox_TS_MEAN=zeros(length(RS_ROI_vox_TS(:,1)),1);
RS_vox_TS_MEAN=zeros(length(RS_vox_TS(:,1)),1);
CORR=zeros(length(RS_ROI_vox_TS(:,1)),length(RS_vox_TS(:,1)));

for i=1:length(RS_ROI_vox_TS(:,1))
    for j=1:length(RS_vox_TS(:,1))
        RS_ROI_vox_TS_MEAN(i,1)=mean(RS_ROI_vox_TS(i,:));
        RS_vox_TS_MEAN(j,1)=mean(RS_vox_TS(j,:));
        % correlation of RS_ROI_vox_TS(i,:) with RS_vox_TS(j,:)
        CORR_NUM=0; %CORR_DEN11=0; CORR_DEN22=0;
        for k=1:length(RS.img(1,1,1,:))
           CORR_NUM1=(RS_ROI_vox_TS(i,k)-RS_ROI_vox_TS_MEAN(i,1))*(RS_vox_TS(j,k)-RS_vox_TS_MEAN(j,1));
           CORR_NUM=CORR_NUM+CORR_NUM1;
           
           %CORR_DEN1=(RS_ROI_vox_TS(i,k)-RS_ROI_vox_TS_MEAN(i,1))^2; CORR_DEN2=(RS_vox_TS(j,k)-RS_vox_TS_MEAN(j,1))^2;
           %CORR_DEN11=CORR_DEN11+CORR_DEN1;
           %CORR_DEN22=CORR_DEN22+CORR_DEN2; %% to find standard deviation manually
        end
        %CORR_DEN=CORR_DEN11*CORR_DEN22;
        %if CORR_DEN~=0
        %    CORR=CORR_NUM/(CORR_DEN11*CORR_DEN22); 
        %end
        CORR_DEN=length(RS.img(1,1,1,:))*StdDev_ROI_TS(i,1)*StdDev_RS_TS(j,1);
        CORR(i,j)=CORR_NUM/CORR_DEN;
    end
end

% Sqrt_CORR=sqrt(CORR); %taking square root of correlation Coefficients
% Sqrt_CORR1=real(Sqrt_CORR); Sqrt_CORR2=-imag(Sqrt_CORR); %seperating into real and imaginary (as sqrt(-1)=imag)
% Sqrt_CC=Sqrt_CORR1+Sqrt_CORR2; %concatenate real and imaginary value 
% CORR_ng=CORR<0; CORR_ps=CORR>=0;
% CORR_neg=CORR_ng.*CORR; CORR_pos=CORR_ps.*CORR;

%%if u wish to save correlation map in .mat format then uncomment next line
% save([feat_loc,'/',ROI_name,'_',num2str(mask_threshold),'/Correlation_result','CORR','-double'])

%% Avg and Max of Correlation Coef.'s (CC) and of Squareroot of CC's of all ROI voxel's and their mapping
RS_MASK_xyzT=RS_MASK_xyz';

 %for avg map
if isempty(find(type_of_map==1,1))==0
    Avg_CC=nanmean(CORR);%mean(CORR,1,'omitnan');
    Avg_CC_MAP=RS_MASK.img;
    for i=1:length(Avg_CC(1,:))
       Avg_CC_MAP(RS_MASK_xyzT(1,i),RS_MASK_xyzT(2,i),RS_MASK_xyzT(3,i))=Avg_CC(1,i);
    end
    
    Avg_CC_img.hdr=RS_MASK.hdr;
    Avg_CC_img.img=Avg_CC_MAP;
    save_nii(Avg_CC_img,[output_dir,'/Avg_CC_map.nii.gz']);
    % registration of map to standard space
    avgmap2std=['sh -c ". ${FSLDIR}etc/fslconf/fsl.sh;${FSLDIR}bin/flirt -in ',output_dir,'/Avg_CC_map.nii.gz',' -ref ','${FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz',' -applyxfm -init ',ROI_dir_name,'/filtered_func2standard.mat',' -out ',output_dir,'/Avg_CC_map_std.nii.gz',' -omat ',output_dir,'/Avg_CC_map_2_std.mat"'];
    system(avgmap2std)
    
    func_mask2std=['sh -c ". ${FSLDIR}etc/fslconf/fsl.sh;${FSLDIR}bin/flirt -in ',feat_loc,'/mask.nii.gz',' -ref ','${FSLDIR}/data/standard/MNI152_T1_1m m_brain.nii.gz',' -applyxfm -init ',output_dir,'/Avg_CC_map_2_std.mat',' -out ',output_dir,'/MASK_Avg_CC_map_std.nii.gz',' -omat ',output_dir,'/MASK_Avg_CC_map_2_std.mat"'];
    system(func_mask2std);
    %load standard space data of subject to construct 4D map
    Std_CC_Avg=load_untouch_nii([output_dir,'/Avg_CC_map_std.nii.gz']);
    Std_CC_Avg_img=Std_CC_Avg.img;
else
    Std_CC_Avg_img=[];
end

% for max map
if isempty(find(type_of_map==2,1))==0
    Max_CC=nanmax(CORR);%max(CORR,[],1,'omitnan');
    Max_CC_MAP=RS_MASK.img;
    for i=1:length(Max_CC(1,:))
       Max_CC_MAP(RS_MASK_xyzT(1,i),RS_MASK_xyzT(2,i),RS_MASK_xyzT(3,i))=Max_CC(1,i);
    end
    
    
    Max_CC_img.hdr=RS_MASK.hdr;
    Max_CC_img.img=Max_CC_MAP;
    save_nii(Max_CC_img,[output_dir,'/Max_CC_map.nii.gz']);
    
    % registration of map to standard space
    maxmap2std=['sh -c ". ${FSLDIR}etc/fslconf/fsl.sh;${FSLDIR}bin/flirt -in ',output_dir,'/Max_CC_map.nii.gz',' -ref ','${FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz',' -applyxfm -init ',ROI_dir_name,'/filtered_func2standard.mat',' -out ',output_dir,'/Max_CC_map_std.nii.gz',' -omat ',output_dir,'/Max_CC_map_2_std.mat"'];
    system(maxmap2std)
    
    func_mask2std=['sh -c ". ${FSLDIR}etc/fslconf/fsl.sh;${FSLDIR}bin/flirt -in ',feat_loc,'/mask.nii.gz',' -ref ','${FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz',' -applyxfm -init ',output_dir,'/Max_CC_map_2_std.mat',' -out ',output_dir,'/MASK_Max_CC_map_std.nii.gz',' -omat ',output_dir,'/MASK_Max_CC_map_2_std.mat"'];
    system(func_mask2std);
    %load standard space data of subject to construct 4D map    
    Std_CC_Max=load_untouch_nii([output_dir,'/Max_CC_map_std.nii.gz']);
    Std_CC_Max_img=Std_CC_Max.img;
else
    Std_CC_Max_img=[];
end

% for avgofmax map
if isempty(find(type_of_map==3,1))==0
    Sort_CORR=sort(CORR,1,'descend');
    AvgofMax_CC=mean(Sort_CORR(1:floor((length(CORR(:,1)))/10),:),1); %Average of maximum ~10% correlation coeficient values
    AvgofMax_CC_MAP=RS_MASK.img;
    for i=1:length(AvgofMax_CC(1,:))    
      AvgofMax_CC_MAP(RS_MASK_xyzT(1,i),RS_MASK_xyzT(2,i),RS_MASK_xyzT(3,i))=AvgofMax_CC(1,i);
    end
    
    
    AvgofMax_CC_img.hdr=RS_MASK.hdr;
    AvgofMax_CC_img.img=AvgofMax_CC_MAP;
    save_nii(AvgofMax_CC_img,[output_dir,'/AvgofMax_CC_map.nii.gz']);
    
    % registration of map to standard space
    avgmap2std=['sh -c ". ${FSLDIR}etc/fslconf/fsl.sh;${FSLDIR}bin/flirt -in ',output_dir,'/AvgofMax_CC_map.nii.gz',' -ref ','${FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz',' -applyxfm -init ',ROI_dir_name,'/filtered_func2standard.mat',' -out ',output_dir,'/AvgofMax_CC_map_std.nii.gz',' -omat ',output_dir,'/AvgofMax_CC_map_2_std.mat"'];
    system(avgmap2std)
    
    func_mask2std=['sh -c ". ${FSLDIR}etc/fslconf/fsl.sh;${FSLDIR}bin/flirt -in ',feat_loc,'/mask.nii.gz',' -ref ','${FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz',' -applyxfm -init ',output_dir,'/AvgofMax_CC_map_2_std.mat',' -out ',output_dir,'/MASK_AvgofMax_CC_map_std.nii.gz',' -omat ',output_dir,'/MASK_AvgofMax_CC_map_2_std.mat"'];
    system(func_mask2std);
    %load standard space data of subject to construct 4D map    
    Std_CC_AvgofMax=load_untouch_nii([output_dir,'/AvgofMax_CC_map_std.nii.gz']);
    Std_CC_AvgofMax_img=Std_CC_AvgofMax.img;
else
    Std_CC_AvgofMax_img=[];
end
else
    disp(['inappropriate ROI or mask threshold value for ',file_loc]);
    disp('function operation is teminated')
    
end


%% 
toc