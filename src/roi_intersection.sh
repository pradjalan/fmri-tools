# This shell script creates an alternate ROI as a binary mask of the voxlels of the input image which lie in the voxels covered in the ROI image. 
# The ROI files are multiplied with the input file. Separate output for positively activated voxels and negatively activated voxels is created.
# In addition to this, the code does not produce any output if there are no voxels in the intersection.

# $1: Input NIfTI file ; $2: Directory Containing ROIs (NIfTI files) ; $3: Image Threshold $4: Output directory for interestion results


for roi in `ls $2 | grep nii`; do

	#Make Directories if they do not exist
	mkdir -p $4

	echo "Doing ROI: $roi"

	#Mask for Positively voxels with positive values. 
	fslmaths $1 -thr $3 -mul $2/$roi  -bin $4/positive_$roi  
	nzm=`fslstats $4/positive_$roi -M`
	
	if [ `echo $nzm == 0 | bc ` == 1 ]; then
		rm $4/positive_$roi
	fi

	#Mask for Positively voxels with negative values. 
	fslmaths $1 -mul -1 -thr $3 -mul $2/$roi -bin $4/negative_$roi 
	nzm=`fslstats $4/negative_$roi -M`
	if [ `echo $nzm == 0 | bc ` == 1 ]; then
		rm $4/negative_$roi
	fi


done