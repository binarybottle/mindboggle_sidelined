#!/usr/bin/python

"""
Transform FreeSurfer surfaces to the original space of a given image.

Copyright 2011 Arno Klein (arno@mindboggle.info) after
               https://gist.github.com/1459725 by Satrajit Ghosh
Apache License, Version 2.0

"""

import os, sys
import numpy as np
import nibabel as nb
from tvtk.api import tvtk
from mayavi import mlab
import nibabel.gifti as gifti

transform_volume = 0
convert_volumes = 0
get_transforms = 1
transform_surfaces = 1
plot_output = 1
plot_output_files = 1

surfaces = ['pial'] #,'smoothwm'] #['inflated','sphere']
hemis = ['lh','rh']

# Inputs
if len(sys.argv) < 3:
    print "Usage: python freesurfer_to_native_space.py <path to FreeSurfer subject> <output directory>"
    exit(-1)
else: 
    subject_path = sys.argv[1]
    output_path = sys.argv[2]

original_volume_mgz = subject_path + '/mri/orig/001.mgz'
conformed_volume_mgz = subject_path + '/mri/brain.mgz'
original_volume_nii = output_path + '/original.nii.gz'
conformed_volume_nii = output_path + '/conformed.nii.gz'
transformed_volume_nii = output_path + '/transformed.nii.gz'

# Reslice a conformed volume in native ("transformed") space:
if transform_volume and os.path.exists(conformed_volume_mgz) and os.path.exists(original_volume_mgz):
    args = ['mri_convert -rl',original_volume_mgz,'-rt nearest',conformed_volume_mgz,transformed_volume_nii]
    print(" ".join(args)); os.system(" ".join(args));  # p = Popen(args); p.close()

# Convert volume files from FreeSurfer's mgz format to nifti format:
if convert_volumes:
    if os.path.exists(original_volume_mgz):
        cmd = 'mri_convert ' + original_volume_mgz + ' ' + original_volume_nii; os.system(cmd)
    if os.path.exists(conformed_volume_mgz):
        cmd = 'mri_convert ' + conformed_volume_mgz + ' ' + conformed_volume_nii; os.system(cmd)

# Load the original and conformed volume files' affine transform matrices:
if get_transforms and os.path.exists(conformed_volume_nii) and os.path.exists(original_volume_nii):
    affine_conformed = nb.load(conformed_volume_nii).get_affine()
    affine_original = nb.load(original_volume_nii).get_affine()

# Create and apply a transform matrix to FreeSurfer's surface meshes:
if transform_surfaces:

    # Create the transform matrix from FreeSurfer's conformed surface to conformed volume:
    M = np.array([[-1,0,0,128],
                  [0,0,1,-128],
                  [0,-1,0,128],
                  [0,0,0,1]],dtype=float)

    print('xfm = np.dot( np.linalg.inv(affine_original), np.dot(affine_conformed, np.linalg.inv(M)))')
    xfm = np.dot( np.linalg.inv(affine_original), np.dot(affine_conformed, np.linalg.inv(M)))

    # Apply the above transform to FreeSurfer's surface meshes:
    for surface in surfaces:
        for hemi in hemis:
            freesurfer_surface = subject_path + '/surf/' + hemi + '.' + surface
            gifti_surface = output_path + '/' + hemi + '.' + surface + '.gii'
            vtk_surface = output_path + '/' + hemi + '.' + surface + '.vtp'
            if os.path.exists(freesurfer_surface):

                # Convert files from FreeSurfer surface format to gifti:
                args = ['mris_convert', freesurfer_surface, gifti_surface]
                print(" ".join(args)); os.system(" ".join(args));  # p = Popen(args); p.close()

                # Load gifti surface:
                surf = gifti.read(gifti_surface)

                # Transform the vertex data array:
                transformed_vertices = np.dot(xfm, np.hstack((surf.darrays[0].data, \
                                       np.ones((surf.darrays[0].data.shape[0],1)))).T)[:3,:].T


                surf=nb.freesurfer.read_geometry(freesurfer_surface)
                xfm = np.dot(affine_conformed, np.linalg.inv(M))
                xfmda0 = np.dot(xfm, np.hstack((surf[0], np.ones((surf[0].shape[0],1)))).T)[:3,:].T
                mesh = tvtk.PolyData(points=xfmda0, polys=surf[1])


                # Create a mesh:
#                mesh = tvtk.PolyData(points=transformed_vertices, polys=surf.darrays[1].data)
                if plot_output and not plot_output_files:
                    mlab.pipeline.surface(mesh)

                # Write gifti surface in original space to vtk file:
                print('w = tvtk.XMLPolyDataWriter(input=mesh, file_name=vtk_surface)')
                w = tvtk.XMLPolyDataWriter(input=mesh, file_name=vtk_surface)
                w.write()

# Plot the brain and surfaces transformed to native space:
if plot_output:
    cdata = nb.load(transformed_volume_nii).get_data()
    #mlab.pipeline.volume(mlab.pipeline.scalar_field(cdata))
    if plot_output_files:
        for surface in surfaces:
            for hemi in hemis:
                 vtk_surface = output_path + '/' + hemi + '.' + surface + '.vtp'
                 mesh_reader = tvtk.XMLPolyDataReader(file_name=vtk_surface)
                 mesh2 = mesh_reader.output
                 mlab.pipeline.surface(mesh2)
    mlab.show()
