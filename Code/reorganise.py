import sys


f = open(sys.argv[1],'r+')
data = f.readlines()
f.close()


additional_conditions = {}
num_add_conds = {}
for line in data:
	if 'alcohol' in line:
		cond = line.split('alcohol')[1].split('(')[0]
		if not cond in additional_conditions:
			additional_conditions[cond] = {}
		if not cond in num_add_conds:
			num_add_conds[cond] = {}
		atts = line.replace('\n').split(';')
		nac = 0
		for att in atts[1:]:
			if not ('Alcohol' in att or 'No Diagnosis' in att):
				nac += 1
				if not att in additional_conditions[cond]:
					additional_conditions[cond][att] = 0
				additional_conditions[cond][att] += 1
		if not nac in num_add_conds[cond]:
			num_add_conds[cond][nac] = 0
		num_add_conds[cond][nac] += 1

print(additional_conditions)
print(num_add_conds)




