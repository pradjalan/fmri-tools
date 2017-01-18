#! /usr/bin/Rscript
require("oro.nifti")

args <- commandArgs(TRUE)

input_file <- args[1]
threshold <- args[2]
atlas_path <- args[3]
out_file <- args[4]

map <-  readNIfTI(input_file, reorient = FALSE)
map_img <- map@.Data

atlas <-  readNIfTI(atlas_path, reorient = FALSE)
atlas_img <- atlas@.Data

active_voxels <- which(map_img>threshold)

out_img <- map_img;
out_img(active_voxels) = atlas_img(active_voxels);

output_nii.hdr=map.hdr;
output_nii.img = out_img;
save_nii(output_nii, out_file);