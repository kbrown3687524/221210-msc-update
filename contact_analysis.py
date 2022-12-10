#!/usr/bin/env python3

from amia import *
import sys
q = MutantAnalysis()
if len(sys.argv) >= 3:
    structures_folder = sys.argv[1]
    model_folder = sys.argv[2]
    cutoff = sys.argv[3]
    output_folder = sys.argv[4]
    q.multiple_residue_contacts(structures_folder, model_folder, cutoff, output_folder)
    q.multi_ligand_contacts_calc(structures_folder, model_folder, output_folder)
elif len(sys.argv) < 3:
	sys.exit('Too few arguments provided')
elif len(sys.argv) > 3:
	sys.exit('Too many arguments provided')
