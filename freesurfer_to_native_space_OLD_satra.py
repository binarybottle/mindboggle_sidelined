import nibabel as nb
import numpy as np
from tvtk.api import tvtk
from mayavi import mlab
import os
#freesurfer surface to volume transform
M = np.array([[-1,0,0,128],
              [0,0,1,-128],
              [0,-1,0,128],
              [0,0,0,1]],dtype=float)
basedir = './KKI2009-01'
surf=nb.freesurfer.read_geometry(os.path.join(basedir, 'surf/lh.pial'))
orig = nb.freesurfer.load(os.path.join(basedir, 'mri/orig.mgz'))
native = nb.freesurfer.load(os.path.join(basedir, 'mri/orig/001.mgz'))
ac = orig.get_affine()
ao = native.get_affine()
xfm = np.dot(ac, np.linalg.inv(M))
xfmda0 = np.dot(xfm, np.hstack((surf[0], np.ones((surf[0].shape[0],1)))).T)[:3,:].T
mesh = tvtk.PolyData(points=xfmda0, polys=surf[1])

native = nb.load('transformed.nii.gz')
native_data = native.get_data()

# Write mesh to vtk file:
vtk_surface = 'lh.pial.native.vtp'
print('w = tvtk.XMLPolyDataWriter(input=mesh, file_name=vtk_surface)')
w = tvtk.XMLPolyDataWriter(input=mesh, file_name=vtk_surface)
w.write()

vol = mlab.pipeline.scalar_field(native_data)
pls = mlab.pipeline.surface(mesh)
plv = mlab.pipeline.volume(vol)
# Thanks to Isaiah Norton:
# The poking is simply telling an actor that you want a transformation applied 
# to your renderer to match it with the physical world.
mt = tvtk.Matrix4x4()
for i in range(4):
    for j in range(4):
        mt.set_element(i,j,ao[i,j])
plv.actors[0].poke_matrix(mt)
plv.update_pipeline()
mlab.show()
