"""
Quick visualisation of a LAMMPS trajectory/data file with OVITO.

GUI route (easiest): open OVITO -> File > Load File -> traj_pH10.lammpstrj (or final_pH10.data).
  * Add "Select type" for type 5 (Ow) + 6 (Hw), then "Delete selected" to hide water.
  * Colour by type: 1 Si (yellow), 2/3/8 O (red), 4 H (white), 7 Na+ (blue), 9 Cl- (green), 10 Ca2+ (pink).
  * Render the last frame (and a frame from the start) to compare ion positions around each particle.

Scripted route (pip install ovito):  python 08_visualize_ovito.py traj_pH10.lammpstrj out.png
"""
import sys
from ovito.io import import_file
from ovito.modifiers import SelectTypeModifier, DeleteSelectedModifier
from ovito.vis import Viewport

infile = sys.argv[1] if len(sys.argv) > 1 else 'traj_pH10.lammpstrj'
outpng = sys.argv[2] if len(sys.argv) > 2 else 'snapshot.png'

pipe = import_file(infile)
pipe.modifiers.append(SelectTypeModifier(property='Particle Type', types={5, 6}))   # water
pipe.modifiers.append(DeleteSelectedModifier())
pipe.add_to_scene()
vp = Viewport(type=Viewport.Type.Perspective, camera_dir=(-1, -1, -1))
vp.zoom_all()
vp.render_image(size=(1600, 900), filename=outpng, frame=pipe.source.num_frames - 1)
print('wrote', outpng)
