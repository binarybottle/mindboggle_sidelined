# Run transform_unconformed_volume_to_standard_space on a directory:

import os

output_path = 'transformed_brainvisa_fundus_volumes_Perrot62'
subjects_path = '/home/arno/Data/Brains/Perrot62_sulci/freesurfer5.1_output_plus_surface_features/'
volumes_path = '/home/arno/Data/Brains/Perrot62_sulci/manually_labeled_brainvisa_fundi/sulci_volumes/'

subjects = os.listdir(subjects_path)
volumes = os.listdir(volumes_path)

for i,volume in enumerate(volumes):
    args = ['python transform_unconformed_volume_to_standard_space.py',
            volumes_path+volume, subjects_path+subjects[i], output_path, volume, 'nearest']
    print(" ".join(args)); os.system(" ".join(args));  # p = Popen(args); p.close()
