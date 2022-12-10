#!/usr/bin/env python3

from amia import *
import sys
p = Mutations()
if len(sys.argv) >= 3:
	mutation_list = sys.argv[1]
	model = sys.argv[2]
	output = sys.argv[3]
	p.process_mutant_list(mutation_list)
	p.multi_mutations(model = sys.argv[2], output = sys.argv[3])
	p.foldx_opti(output)
elif len(sys.argv) < 3:
	sys.exit('Too few arguments provided')
elif len(sys.argv) > 3:
	sys.exit('Too many arguments provided')
