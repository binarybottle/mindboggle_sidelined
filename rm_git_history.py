rm_list = ['doc/_build',
'doc/api/generated',
'doc/bibtex',
'doc/devel',
'doc/images',
'doc/_static/logos',
'docs/_build',
'dti',
'eliezer*'
'extract/fundi_vtk',
'extract/medial_surface',
'extract/medial_surfaces',
'extract/pits',
'extract/spectral',
'extract/surface_depth',
'evaluate/x',
'feature_extraction/basins_pits_im_2011',
'feature_extraction/crests',
'feature_extraction/fundi/li_2010',
'feature_extraction/pits',
'feature_extraction/sift',
'feature_extraction/skeleton',
'feature_extraction/tests',
'feature_matching',
'Makefile*',
'Manifest.in',
'measure/surface_measures',
'measure/surface_travel_depth',
'mindboggle/propagate/realignment_test',
'mindboggle/utils/alternatives.py',
'mindboggle/utils/cp_freesurfer_files.py',
'mindboggle/utils/debug_nipype.py',
'mindboggle/utils/samples',
'pipeline'
'publication',
'propagate/.idea',
'propagate/realignment_test',
'qtcreator-build',
'register',
'registration',
'segment',
'surface_travel_depth',
'tools',
'utilities',
'utils/cp_freesurfer_files.py',
'utils/register/ANTS_Multimodality_notes.txt',
'x',
'.idea',
'*.png',
'*.pyc',
'*~']

import os

#for rm_file in rm_list:
#    os.system("git filter-branch -f --index-filter 'git rm -r --cached --ignore-unmatch " + rm_file + "' --prune-empty --tag-name-filter cat -- --all")
os.system("rm -r .git/refs/original")
os.system("git push origin master --force")
os.system("git reflog expire --expire=now --all")
os.system("git gc --prune=now")
os.system("git gc --aggressive --prune=now")

