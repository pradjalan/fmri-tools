from lxml import etree
import nibabel as nib
import numpy as np 
import os, scipy.io
import copy


atlas_name = 'IITD-CellType_GyrusLobe-TalairachCerebellar'

standard_path = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain.nii.gz'

atlas_path = '/usr/share/fsl/5.0/data/atlases/Talairach/Talairach-labels-2mm.nii.gz'
atlas_xml_path = '/usr/share/fsl/5.0/data/atlases/Talairach.xml'

cerebellar_atlas_path = '/usr/share/fsl/5.0/data/atlases/Cerebellum/Cerebellum-MNIflirt-maxprob-thr50-2mm.nii.gz'
cerebellar_xml_path = '/usr/share/fsl/5.0/data/atlases/Cerebellum_MNIflirt.xml' 

output_path = '/mnt/project1/home1/cs5120287/my_atlases/'+atlas_name+'-labels-2mm.nii.gz'
xml_output = '/mnt/project1/home1/cs5120287/my_atlases/'+atlas_name+'.xml'
csv_directory = '/mnt/project1/home1/cs5120287/my_atlases/'



new_roi_indices = {} #ROI_NAME -> NEW_INDEX
index_dictionary = {} # NEW_INDEX -> OLD_INDICES
old_roi_indices = {} # OLD_INDEX -> NEW_INDEX
cur_index = -1



atlas_tree = etree.parse(atlas_xml_path)
sample_child = etree.tostring(atlas_tree.find('data').getchildren()[0])


#Add ROIs from original Atlas (Talairach)
for child in atlas_tree.find('data').getchildren():
	child_index = int(child.get('index'))
	roi_data = child.text.split('.')
	hemisphere = roi_data[0]
	lobe = roi_data[1]
	gyrus = roi_data[2]
	tissue = roi_data[3]
	cell_type = roi_data[4]

	if 'Cerebellum' in hemisphere:
		continue

	
	roi_name = hemisphere+'.'

	
	

	if cell_type!='*':
		roi_name += '*.*.*.'+cell_type
	else:
		roi_name += lobe + '.' + gyrus + '.*.*'


	
	if not roi_name in new_roi_indices:
		cur_index += 1
		new_roi_indices[roi_name] = cur_index
		print(str(cur_index)+','+roi_name)
	if not new_roi_indices[roi_name] in index_dictionary:
		index_dictionary[new_roi_indices[roi_name]] = []
	index_dictionary[new_roi_indices[roi_name]].append(child_index)
	old_roi_indices[child_index] = new_roi_indices[roi_name]





atlas_tree = etree.parse(cerebellar_xml_path)
sample_child = etree.tostring(atlas_tree.find('data').getchildren()[0])

index_dictionary_cerebellar = {} # NEW_INDEX -> OLD_INDICES
old_roi_indices_cerebellar = {} # OLD_INDEX -> NEW_INDEX

# #Add ROIs from Cerebellar Atlas
for child in atlas_tree.find('data').getchildren():
	child_index = int(child.get('index'))+1
	roi_data = child.text.split(' ')
	hemisphere = roi_data[0] + ' Cerebellum'
	description = ' '.join(roi_data[1:])

	
	roi_name = hemisphere+'.*.'+description+'.*.*'
	


	if not roi_name in new_roi_indices:
		cur_index += 1
		new_roi_indices[roi_name] = cur_index
		print(str(cur_index)+','+roi_name)
	if not new_roi_indices[roi_name] in index_dictionary_cerebellar:
		index_dictionary_cerebellar[new_roi_indices[roi_name]] = []
	index_dictionary_cerebellar[new_roi_indices[roi_name]].append(child_index)
	old_roi_indices_cerebellar[child_index] = new_roi_indices[roi_name]




atlas = nib.load(atlas_path)
cerebellar_atlas = nib.load(cerebellar_atlas_path)

atlas_data = atlas.get_data()
cerebellar_atlas_data = cerebellar_atlas.get_data()

# my_atlas_data = np.empty(atlas_data.shape,dtype=int)
my_atlas_data = copy.deepcopy(atlas_data)
my_atlas_data *= 0

for i in old_roi_indices:
	my_atlas_data[np.where(atlas_data==i)] = old_roi_indices[i]
for i in old_roi_indices_cerebellar:
	my_atlas_data[np.where(cerebellar_atlas_data==i)] = old_roi_indices_cerebellar[i]



my_atlas = nib.Nifti1Image(my_atlas_data, np.eye(4), atlas.header)
my_atlas.header.set_dim_info(2,2,2)
nib.save(my_atlas,output_path)
scipy.io.savemat(csv_directory+'/'+atlas_name+'.mat', mdict={'my_atlas_data':my_atlas_data})

#Change Header of the new atlas to be the same as the original atlas.


# os.system('fslchfiletype NIFTI_PAIR '+standard_path+' '+csv_directory+'/temp')
# os.system('fslchfiletype NIFTI_PAIR '+output_path+' '+csv_directory+'/temp_out')
# os.system('cp '+csv_directory+'/temp.hdr ' + csv_directory+'temp_out.hdr')
# os.system('fslchfiletype NIFTI_GZ '+csv_directory+'/temp_out ' + output_path)




new_atlas_tree = etree.fromstring(etree.tostring(atlas_tree))
for c in new_atlas_tree.find('data').getchildren():
	c.getparent().remove(c)
f = open(csv_directory+atlas_name+'_regionIndices.csv', 'w')
for l in new_roi_indices:
	f.write(str(new_roi_indices[l]) + ',' + l + '\n')
	this_element = etree.fromstring(sample_child)
	this_element.set('index',str(new_roi_indices[l]))
	this_element.text = l
	new_atlas_tree.find('data').append(this_element)

f = open(xml_output,'w')
f.write(etree.tostring(new_atlas_tree, pretty_print=True))
f.close()


f= open(csv_directory+atlas_name+'_variables.txt','w')
print >> f,('new_roi_indices='+str(new_roi_indices))
print >> f,('old_roi_indices='+str(old_roi_indices))
print >> f,('index_dictionary='+str(index_dictionary))
f.close()



f = open(csv_directory+atlas_name+'_oldIndices.csv','w')
for i in sorted(old_roi_indices.keys()):
		f.write(str(i) + ',' + str(old_roi_indices[i]) + '\n')
f.close()

f = open(csv_directory+atlas_name+'_oldIndices.csv','a')
for i in sorted(old_roi_indices_cerebellar.keys()):
		f.write(str(i) + ',' + str(old_roi_indices_cerebellar[i]) + '\n')
f.close()


