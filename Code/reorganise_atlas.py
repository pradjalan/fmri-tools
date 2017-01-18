from lxml import etree
import nibabel as nib
import numpy as np 


atlas_path = '/usr/share/fsl/5.0/data/atlases/Talairach/Talairach-labels-2mm.nii.gz'
atlas_xml_path = '/usr/share/fsl/5.0/data/atlases/Talairach.xml'
output_path = '/mnt/project1/home1/cs5120287/my_atlases/Modified-Coarse-Talairach-labels-2mm.nii.gz'
xml_output = '/mnt/project1/home1/cs5120287/my_atlases/Modified-Coarse-Talairach.xml'
csv_directory = '/mnt/project1/home1/cs5120287/my_atlases/'


new_roi_indices = {} #ROI_NAME -> NEW_INDEX
index_dictionary = {} # NEW_INDEX -> OLD_INDICES
old_roi_indices = {} # OLD_INDEX -> NEW_INDEX
cur_index = -1



atlas_tree = etree.parse(atlas_xml_path)
sample_child = etree.tostring(atlas_tree.find('data').getchildren()[0])

for child in atlas_tree.find('data').getchildren():
	child_index = child.get('index')
	roi_data = child.text.split('.')
	roi_name = roi_data[0]+'.'+roi_data[2]#+'.'+roi_data[4]
	# if '*.*' in roi_name:
	if '*' in roi_name:
		roi_name = roi_data[0]+'.'+roi_data[1]

	if not roi_name in new_roi_indices:
		cur_index += 1
		new_roi_indices[roi_name] = cur_index
	if not cur_index in index_dictionary:
		index_dictionary[cur_index] = []
	index_dictionary[cur_index].append(child_index)
	old_roi_indices[child_index] = cur_index


	


atlas = nib.load(atlas_path)
my_atlas_data = atlas.get_data()
for i in old_roi_indices:
	my_atlas_data[np.where(my_atlas_data==i)] = old_roi_indices[i]
my_atlas = nib.Nifti1Image(my_atlas_data, np.eye(4))
nib.save(my_atlas,output_path)


new_atlas_tree = etree.fromstring(etree.tostring(atlas_tree))
for c in new_atlas_tree.find('data').getchildren():
	c.getparent().remove(c)
f = open(csv_directory+'/Modified-Coarse-Talairach_regionIndices.csv', 'w')
for l in new_roi_indices:
	f.write(str(new_roi_indices[l]) + ',' + l + '\n')
	this_element = etree.fromstring(sample_child)
	this_element.set('index',str(new_roi_indices[l]))
	this_element.text = l
	new_atlas_tree.find('data').append(this_element)

f = open(xml_output,'w')
f.write(etree.tostring(new_atlas_tree, pretty_print=True))
f.close()



