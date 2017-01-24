function clusterize(input_file, lower_threshold, upper_threshold, atlas_path, out_file)
    %Assuming that the input file is already in the same domain as the atlas file 
    
    map = load_untouch_nii(input_file);
    map_img = map.img;
    
    atlas = load_untouch_nii(atlas_path);
    atlas_img = atlas.img;
    
    active_voxels = find(map_img>lower_threshold & map_img<upper_threshold);
    
    map_img = zeros(size(atlas),'int16');
    map_img(active_voxels) = atlas_img(active_voxels);
    
    
    map.img = map_img; 
    save_untouch_nii(map, out_file);
end
