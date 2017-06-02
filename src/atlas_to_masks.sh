# $1: input file, $2: min intensity $3: max intensity, $4: output dir; $5: Atlas prefix

mkdir -p $4
for i in `seq $2 $3`; do
	echo $i
	n="$i_$5.nii.gz"
	echo $n
	fslmaths $1 -thr $i -uthr $i -bin $4/$n
done