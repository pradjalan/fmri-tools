function fdr(p_file,t_file,mask_file, out_file_path)

% This Code is a part of the fmri-tools utilities

% fdr_correction(p_fileLoc,q_value,out_dir)
% Creates a fdr corrected file with -log(q-values)(to the base 10) using the <input file> which contains -log(p-values) (to the base 10) for the voxels 
% (negative values are permitted and signs are copied onto output q-values). 
% If no mask or threshold is specified then only the non-zero voxels in the input are considered for testing, others are ignored. 
% In case mast file is specified, voxels with non-zero values in the mask file are used for testing. 
% In case the <input file> contains multiple volumes, all the volumes are assumed to be p-values and are converted to q-values if --vol option is not given. 
% If --vol option is used then only the specified <volume number> is converted from p-values to q-values, the rest of the volumes are copied to the output.
%%FDR correction

% If t_file is also mentioned, this code will apply the mask of p-values that pass the FDR correction procedure

fsldir = getenv('FSLDIR');

p_map_nii=load_untouch_nii(p_file);


p_map_img=p_map_nii.img; %seperating img file from nii file


if length(out_file_path)==0
    [out_dir,name,extension]=fileparts(p_file);
    out_file_name=[name extension]; % specify output file name
else
    [out_dir,name,extension]=fileparts(out_file_path);
    out_file_name=[name extension]; % specify output file name
end

if length(mask_file)==0
    mask_file = [fsldir,'/data/standard/MNI152_T1_2mm_brain_mask.nii.gz'];    
end    

std_brain_mask=load_untouch_nii(mask_file);
std_mask_img=double(std_brain_mask.img); %seperating img file from standard nii file



mask_indices=find(std_mask_img); %to find indices of brain mask

p_values=p_map_img(mask_indices); %extracting p_values of brain only (i.e. excluding outside brain region)
sign_mask = sign(p_values);



[sorted_p_values,sorted_p_indices] = sort(abs(p_values),'descend');

q_star_values = zeros(size(sorted_p_values));
rejected = zeros(size(sorted_p_values));

LENGTH=length(sorted_p_values);

q_value = 0.01;

%plot_data
% qis = zeros(length(sorted_p_values);
xaxis = (1:LENGTH)/LENGTH;
qis = q_value*xaxis;


%note that the q_values are to be stored as -log(q) and p values are
%assumed to be stored as -log(p)

%FDR Correction and BH FDR control
for i=1:length(sorted_p_values)
     %Finding q-value
     q_own = sorted_p_values(i) + log10(i) - log10(LENGTH);
     q_star_values(i) = q_own;     

     %fill plot data
%      qis(i) = q_value*(i/LENGTH);
     
     % BH Control Procedure
     if sorted_p_values(i) >= -log10(i) + log10(LENGTH) - log10(q_value)
         rejected(i) = 1;
         continue
     end
end

%FDR Adjustment
q_min = -inf;
for i=abs(-length(sorted_p_values):-1)
    if q_star_values(i) < q_min
        q_star_values(i) = q_min;
    else
        q_min = q_star_values(i);
    end
end




%Consider Rejected Hypothesis
rejected_indices = find(rejected);

p_val_mapped_on_brain=zeros(size(p_values));
p_val_mapped_on_brain(:) = 0;
p_val_mapped_on_brain(sorted_p_indices(rejected_indices)) = sorted_p_values(rejected_indices);

q_star(sorted_p_indices) = q_star_values;
disp(['Voxels with Non-Zero Correlation: ' num2str(100*length(find(abs(q_star)>2))/length(mask_indices)) '%']);


%Create Mask of p-values with rejected hypothesis i.e. a mask of selected p-values
selected_p_values = false(size(mask_indices));
selected_p_values(sorted_p_indices(rejected_indices)) = 1;

disp(['Voxels with Non-Zero Correlation (rejected hypothesis): ' num2str(100*length(find(selected_p_values))/length(mask_indices)) '%']);


%Preserve the negative sign for negative p-values
p_val_mapped_on_brain = sign_mask.*p_val_mapped_on_brain;
q_star = sign_mask.*q_star.';

% Write FDR Corrected p-Values and q-Values to output files
std_mask_img(mask_indices)=p_val_mapped_on_brain;
fdr_corr_nii.img=std_mask_img; fdr_corr_nii.hdr=p_map_nii.hdr;
save_nii(fdr_corr_nii,[out_dir '/fdr_p_values_' out_file_name]);

std_mask_img(mask_indices)=q_star;
q_fdr_nii.img = std_mask_img; q_fdr_nii.hdr = p_map_nii.hdr;
save_nii(q_fdr_nii,[out_dir '/q_values_' out_file_name]);

std_mask_img(mask_indices) = selected_p_values;
fdr_mask_nii.img = std_mask_img; fdr_mask_nii.hdr = p_map_nii.hdr;
save_nii(fdr_mask_nii,[out_dir '/fdr_mask_' out_file_name]);


%plot figures
figure('visible','off');
hold on;
plot(xaxis,qis,'r');
plot(xaxis,(10^(-sorted_p_values)),'b--o');
legend('q*(i/N)','p_values');
xlabel('i/N');
ylabel('p_value');
title(['FDR Correction: ' out_dir]);

saveas(gcf,[out_dir '/fdr_plot_' name '.png']);



%Apply Mask to the t-values
if length(t_file)~=0
    ApplyMask = ['sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/fslmaths ',t_file,' -mul ',out_dir,'/fdr_mask_',out_file_name,'  ',out_dir,'/fdr_t_values_',out_file_name,'  "'];
    sam = system(ApplyMask);
end

disp(['Saved FDR corrected files in: ',out_dir])
