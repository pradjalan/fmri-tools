# $1:directory list, $2: TR identification code=645|1400|CAP, $3:TR value(sec)=0.645|1.4|2.5, $4:outpath (base), $5:from, $6 to

count=0
source_path="/mnt/project1/rawData/fMRI/incoming/BrainScape_fBIRN/NKI-RS-Lite/COINS/Data/dicom2nii/"
for dr in `cat $1 | grep /REST_$2 `; do
	if [ $count -le $6 ]; then
	if [ $5 -le $count ]; then	
	
	if [ ! -d $4/$dr/ ]; then
            mkdir -p $4/$dr/
    fi
    echo "Merging: $dr"
	fslmerge  -tr $4/$dr/rest.nii $source_path/$dr/*.nii $3
	echo "Scans Merged: $count.."

	fi
	fi
	count=$[$count +1]
done