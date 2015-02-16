# convert LINE segment in VTK into SCALAR LUT, for any vertex shown in LINE, it has a value greater than -1 in a scalar LUT called fundi 
# For Forrest's fundi and Gang Li's fundi

def load_vtk(File):
    """Load the VTK file to get faces, points, lines and number of points
    """
    from mindboggle.utils.io_vtk import read_vtk
    faces, lines, u3, points, npoints, u6, u7, u8 = read_vtk(File)
    return faces, points, lines, npoints

def line2scalar(lines, npoints):
    """Convert LINES in VTK into an SCALAR LUT
    """
    import itertools
    fundi = set(list(itertools.chain(*lines)))
    scalar = [1 if vertex in fundi else -1 for vertex in range(1, npoints+1)]
    return scalar

def rewrite_vtk(Faces, points, scalar, output_vtk):
    """write scalar-represented fundi file
    """
    from mindboggle.utils.io_vtk import write_vtk
    write_vtk(output_vtk, points, faces=Faces, scalars=[scalar], scalar_names=['fundi'])

def loop_thru_forrest(data_path):
    """Loop thru all subjects under data_path

    This will OVERWRITE old fundi file! 

    """
    import os.path
    import os
    for hemi_subj in os.listdir(data_path):
        if "lh" in hemi_subj:
            line_fundi_vtk = os.path.join(data_path, hemi_subj, "lh.pial.fundi.vtk")
            scalar_fundi_vtk = os.path.join("/data/data/Mindboggle_MRI/MB101/results/Forrest_Fundi.scalars", hemi_subj, "lh.pial.fundi.vtk")
        elif "rh" in hemi_subj:
            line_fundi_vtk = os.path.join(data_path, hemi_subj, "rh.pial.fundi.vtk")
            scalar_fundi_vtk = os.path.join("/data/data/Mindboggle_MRI/MB101/results/Forrest_Fundi.scalars", hemi_subj, "rh.pial.fundi.vtk")
        else:
            print "wrong subject or wrong hemisphere"

        print "=> ", line_fundi_vtk
        print "<= ", scalar_fundi_vtk

        faces, points, lines, npoints = load_vtk(line_fundi_vtk)
        scalar = line2scalar(lines, npoints)

        rewrite_vtk(faces, points, scalar, scalar_fundi_vtk)

def loop_thru_gangli(data_path):
    """Loop thru all subjects under data_path
    For Gang Li's fundi 

    There are two kinds of fundi from Gang Li's code, *.SulciPruneCombine.vtk and *.SulciPruneCombineReconnect.vtk 
    If *.SulciPruneCombineReconnect.vtk exists, do not use *.SulciPruneCombine.vtk.
    Otherwise, use it. 

    This function is not usable as Gang Li's VTK has only fundus vertexes. 

    """

    import os.path
    import os
    for hemi_subj in os.listdir(data_path):
        if "lh" in hemi_subj:
            hemi = "lh"
        elif "rh" in hemi_subj:
            hemi = "rh"
        else:
            print "wrong subject or wrong hemisphere"

        line_fundi_vtk = os.path.join(data_path, hemi_subj, hemi+".pial.SulciPruneCombineReconnect.vtk")

        if os.path.exists(line_fundi_vtk):
            pass
        else:
            line_fundi_vtk = os.path.join(data_path, hemi_subj, hemi+".pial.SulciPruneCombine.vtk")
            if os.path.exists(line_fundi_vtk):
                pass
            else:
                print ("cannot locate fundi file")

        scalar_fundi_vtk = os.path.join(data_path, hemi_subj, hemi+".pial.fundi.vtk")

        print "=> ", line_fundi_vtk
        print "<= ", scalar_fundi_vtk

        faces, points, lines, npoints = load_vtk(line_fundi_vtk)
        scalar = line2scalar(lines, npoints)
        rewrite_vtk(faces, points, scalar, scalar_fundi_vtk)

if __name__ == "__main__":
#    faces, points, lines, npoints = load_vtk("/data/data/Mindboggle_MRI/MB101/results/Forrest_Fundi.lines/_hemi_lh_subject_Twins-2-1/lh.pial.fundi.vtk")
#    scalar = line2scalar(lines, npoints)
#    output_vtk = "test.vtk"
#    rewrite_vtk(faces, points, scalar, output_vtk)

    loop_thru_forrest("/data/data/Mindboggle_MRI/MB101/results/Forrest_Fundi.lines")
#    loop_thru_gangli("/data/data/Mindboggle_MRI/MB101/results/GangLiFundi")
