#!/usr/bin/env python3

from amia import *
import sys
m= FoldXAnalysis()
if len(sys.argv) >= 3:
    structures_folder = sys.argv[1]
    model_path = sys.argv[2]
    output_folder = sys.argv[3]
    foldx_exe = sys.argv[4]
    m.foldx_stability(structures_folder , model_path, output_folder, foldx_exe)
elif len(sys.argv) < 3:
	sys.exit('Too few arguments provided')
elif len(sys.argv) > 4:
	sys.exit('Too many arguments provided')
