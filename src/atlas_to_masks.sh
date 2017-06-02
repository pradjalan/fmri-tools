# $1: input file, $2: min intensity $3: max intensity, $4: output dir


for i in `seq $2 $3`; do
	echo $i
	fslmaths $1 -thr $i -uthr $i -bin $4/$i_$1
done