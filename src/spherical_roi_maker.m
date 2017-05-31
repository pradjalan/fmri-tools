function spherical_roi_maker(config_file,out_dir,mask_file)
%config file contains coordinates, radii and names


fsldir = getenv('FSLDIR');

if length(mask_file)==0
    mask_file = [fsldir,'/data/standard/MNI152_T1_2mm_brain_mask.nii.gz'];
end

if ~exist(out_dir,'dir')
    system(['mkdir -p ' out_dir]);
end
standard_header=load_untouch_header_only(mask_file);
standard_file = load_untouch_nii(mask_file);
standard_image = standard_file.img;

fid = fopen(config_file);


xs = str2num(fgetl(fid));
ys = str2num(fgetl(fid));
zs = str2num(fgetl(fid));
radii = str2num(fgetl(fid));
names = strsplit(fgetl(fid));
use_mm = str2num(fgetl(fid));

xyzs = [xs; ys; zs];
% xyzs = xyzs.';
disp(size(xyzs));

if use_mm
    [~,voxel_coordinates] = map_coords(xyzs,mask_file);
else
    voxel_coordinates = xyzs;
end
voxel_coordinates = floor(voxel_coordinates);


fclose(fid);

for cur=1:length(names)
 roi = char(names(cur));
 xyz = voxel_coordinates(:,cur);
 x = xyz(1);
 y = xyz(2);
 z = xyz(3);
 radius = radii(cur);

    disp(['ROI: ' roi ';  Coordinates: ' num2str(x) ',' num2str(y) ',' num2str(z) ';  ' 'Radius: ' num2str(radius)]);
 
  % Dimensions of standard mask file
%      roi_mask_image = zeros(size(standard_header.dime.dim(2:4)));
    roi_mask_image = standard_image;
    roi_mask_image(find(roi_mask_image)) = 0;
 
    for i=(x-radius):(x+radius)
      for j=(y-radius):(y+radius)
          for k=(z-radius):(z+radius)
            if abs(i-x)+abs(j-y)+abs(k-z) <= radius
                  roi_mask_image(i,j,k) = 1;
            end
          end
      end
     end
 
 roi_mask.hdr = standard_header;
 roi_mask.img = roi_mask_image;
 save_nii(roi_mask, [out_dir '/' roi '_' num2str(radius) '.nii.gz']);
 
end
     
% roi_mask_image = zeros(standard_header.dime.dim(2:4));
% standard_header.dime.dim(2:4);
% sz = size(roi_mask_image);
% 
% roi_mask_image(1,1,1) = 111;
% roi_mask_image(sz(1),1,1) = 211;
% roi_mask_image(1,sz(2),1) = 121;
% roi_mask_image(1,1,sz(3)) = 112;
% roi_mask_image(sz(1),sz(2),1) = 221;
% roi_mask_image(sz(1),1,sz(3)) = 212;
% roi_mask_image(1,sz(2),sz(3)) = 122;
% roi_mask_image(sz(1),sz(2),sz(3)) = 222;
% 
% 
% roi_mask.hdr = standard_header;
% roi_mask.img = roi_mask_image;
% save_nii(roi_mask, [out_dir '/' 'a_1' '.nii.gz']);


end