import sys, os

search_dir_file = sys.argv[1]
subjects_file = open(sys.argv[2],'r+')
out_file = sys.argv[3]




for line in subjects_file.readlines():
	line = line.replace('\n','')
	os.system('cat '+search_dir_file+' | grep '+ line +' >> '+out_file)

subjects_file.close()

