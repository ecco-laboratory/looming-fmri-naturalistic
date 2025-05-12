# %%
# imports
# you must activate the conda env /home/data/eccolab/CondaEnvs/pymeshlab to run this script
import os
import pymeshlab
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

base_mesh_dir = '/archival/projects/SPLaT/data/mesh'
# %%
# Argle parser
parser = ArgumentParser(
    prog='combine_downsample_recon_mesh',
    description='Use pymeshlab to combine and downsample a bilateral freesurfer cortical mesh for 3D printing.',
    formatter_class=ArgumentDefaultsHelpFormatter
)
parser.add_argument(
    '--subj_num',
    '-s',
    type=int,
    help="Subject number to process. Only the number, no leading 0s."
)

args = vars(parser.parse_args())

subj_num = args['subj_num']

# %%
# Do the dirty deed
this_mesh_dir = os.path.join(base_mesh_dir, 'sub-{:04d}'.format(subj_num), 'surf')
mesh_bl = pymeshlab.MeshSet()
mesh_bl.load_new_mesh(os.path.join(this_mesh_dir, 'lh.stl'))
mesh_bl.load_new_mesh(os.path.join(this_mesh_dir, 'rh.stl'))
# the default working directory within a called python script is the current terminal dir
# so absolute path the filter script so you can call this runner script from wherever
mesh_bl.load_filter_script('/home/data/eccolab/SPLaT_fMRI/code/python/meshlab_simplify_brain.mlx')
mesh_bl.apply_filter_script()
mesh_bl.save_current_mesh(
    file_name = os.path.join(this_mesh_dir, 'bl.stl'),
    save_face_color = False
    )
