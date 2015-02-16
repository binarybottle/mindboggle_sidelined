# Script for running Mindboggle on the Mindboggle-101 set:
import os
from mindboggle.utils.io_table import read_columns
from mindboggle.utils.utils import execute

atlas_list_file = '/brains/Mindboggle101/code/mindboggle101_atlases.txt'
atlas_list = read_columns(atlas_list_file, 1)[0]

c = ' '.join(['python mindboggler.py', '--sulci', '--fundi', '--vertices',
              '--spectra 10', '--moments 10', '--visual hier', '--surface_labels manual'])
 
for atlas in atlas_list:
    if 'MMRR-21-' in atlas:
        cmd = c + ' ' + atlas
        execute(cmd)

for atlas in atlas_list:
    if 'NKI-RS-' in atlas:
        cmd = c + ' ' + atlas
        execute(cmd)

for atlas in atlas_list:
    if 'NKI-TRT-' in atlas:
        cmd = c + ' ' + atlas
        execute(cmd)

for atlas in atlas_list:
    if 'OASIS-TRT-20' in atlas:
        cmd = c + ' ' + atlas
        execute(cmd)

for atlas in atlas_list:
    if 'HLN-' in atlas or 'Twins-' in atlas or 'Colin' in atlas or 'After' in atlas or '3T7T' in atlas:
        cmd = c + ' ' + atlas
        execute(cmd)

