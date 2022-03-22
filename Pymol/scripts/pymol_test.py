import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI

import sys, time, os
import pymol

pymol.finish_launching()

##
# Read User Input
spath = os.path.abspath(sys.argv[1])
sname = spath.split('/')[-1].split('.')[0]

# Load Structures

pymol.cmd.load(spath, sname)
pymol.cmd.disable("all")
pymol.cmd.enable(sname)
pymol.cmd.png("my_image.png")

# Get out!
pymol.cmd.quit()