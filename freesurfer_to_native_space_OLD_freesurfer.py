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
import nibabel.gifti as gifti
from tvtk.api import tvtk

transform_volume = 1
convert_volumes = 1
get_transforms = 1
transform_surfaces = 1
plot_volume_surfaces = 1

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
original_volume_mhd = output_path + '/original.mhd'
conformed_volume_nii = output_path + '/conformed.nii.gz'
transformed_volume_mgz = output_path + '/transformed.mgz'
transformed_volume_nii = output_path + '/transformed.nii.gz'
transformed_volume_mhd = output_path + '/transformed.mhd'
transform = output_path + '/transform.dat'
vtk_surface_left = output_path + '/lh.pial.vtp'
vtk_surface_right = output_path + '/rh.pial.vtp'

# Reslice a conformed volume in native ("transformed") space:
if transform_volume and os.path.exists(conformed_volume_mgz) and os.path.exists(original_volume_mgz):
    """
    args = ['tkregister2 --mov',conformed_volume_mgz,'--targ',original_volume_mgz,'--regheader --noedit --reg',transform]
    print(" ".join(args)); os.system(" ".join(args));  # p = Popen(args); p.close()
    args = ['mri_vol2vol --mov',conformed_volume_mgz,'--targ',original_volume_mgz,'--reg',transform,'--o',transformed_volume_mgz]
    print(" ".join(args)); os.system(" ".join(args));  # p = Popen(args); p.close()
    """
    args = ['mri_convert -rl',original_volume_mgz,'-rt nearest',conformed_volume_mgz,transformed_volume_mgz]
    print(" ".join(args)); os.system(" ".join(args));  # p = Popen(args); p.close()

# Convert volume files from FreeSurfer's mgz format to nifti format:
if convert_volumes:
    if os.path.exists(original_volume_mgz):
        cmd = 'mri_convert ' + original_volume_mgz + ' ' + original_volume_nii; os.system(cmd)
    if os.path.exists(conformed_volume_mgz):
        cmd = 'mri_convert ' + conformed_volume_mgz + ' ' + conformed_volume_nii; os.system(cmd)
    if os.path.exists(transformed_volume_mgz):
        cmd = 'mri_convert ' + transformed_volume_mgz + ' ' + transformed_volume_nii; os.system(cmd)

# Load the original and conformed volume files' affine transform matrices:
if get_transforms and os.path.exists(conformed_volume_nii) and os.path.exists(original_volume_nii):
    affine_conformed = nb.load(conformed_volume_nii).get_affine()
    affine_original = nb.load(original_volume_nii).get_affine()

# Create and apply a transform matrix to FreeSurfer's surface meshes:
if transform_surfaces:

    # Create the transform matrix from FreeSurfer's conformed volume to conformed surface:
    M1 = np.array([[-1,0,0,128],
                   [0,0,1,-128],
                   [0,-1,0,128],
                   [0,0,0,1]],dtype=float)

    # Create the transform matrix to translate from FreeSurfer's original volume to vtk original volume:
    M2 = np.array([[1,0,0,218],
                   [0,1,0,0],
                   [0,0,1,0],
                   [0,0,0,1]],dtype=float)

    # Create a transform matrix to go from Freesurfer's conformed surface to original space in VTK coordinates.
    print('xfm = np.dot( M2, np.dot( np.linalg.inv(affine_original), np.dot(affine_conformed, np.linalg.inv(M1))))')
    xfm = np.dot( M2, np.dot( np.linalg.inv(affine_original), np.dot(affine_conformed, np.linalg.inv(M1))))

    # Apply the above transform to FreeSurfer's surface meshes:
    surfaces = ['pial'] #,'smoothwm'] #['inflated','sphere']
    hemis = ['lh','rh']
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

                # Write gifti surface in original space to vtk file:
                mesh = tvtk.PolyData(points=transformed_vertices, polys=surf.darrays[1].data)
                print('w = tvtk.XMLPolyDataWriter(input=mesh, file_name=vtk_surface)')
                w = tvtk.XMLPolyDataWriter(input=mesh, file_name=vtk_surface)
                w.write()

# Plot the original volume and transformed surfaces
if plot_volume_surfaces:

    # Convert original volume file from nifti format to mhd format:
    volume_to_plot_nii = transformed_volume_nii # original_volume_nii
    volume_to_plot_mhd = transformed_volume_mhd # original_volume_mhd

    cmd = 'c3d ' + volume_to_plot_nii + ' -o ' + volume_to_plot_mhd; os.system(cmd)

    args = ['mayavi2 -d', volume_to_plot_mhd, '-m Volume -d', 
                          vtk_surface_left, '-m Surface -d', vtk_surface_right, '-m Surface &']
    print(" ".join(args)); os.system(" ".join(args));  # p = Popen(args); p.close()
