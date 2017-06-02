# $1: input file, $2: min intensity $3: max intensity, $4: output dir; $5: Atlas prefix


for i in `seq $2 $3`; do
	echo $i
	mkdir -p $4
	fslmaths $1 -thr $i -uthr $i -bin $4/$i_$5.nii.gz
done