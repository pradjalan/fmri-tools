import sys
sf = open(sys.argv[1],'r+')
subjects = sf.read().split('\n')
sf.close()

ref = open(sys.argv[2],'r+')

outf = open(sys.argv[3],'w+')
for line in ref:
	sub = line.split('/')[5]
	if sub in subjects:
		outf.write(line)

ref.close()
outf.close()



