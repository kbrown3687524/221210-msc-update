# Project Title: "Automated computational workflow to prioritize potential resistance variants identified in HIV Integrase Subtype C and CRF02_AG" 
# This script is developed for the fufuillment for Masters at the South African National Bioinformatics Institute at the University of the Western Cape
# The project is funded by the Poliomyelitis Research Foundation and the UWC Ada & Bertie Levenstein Bursary Programme
# Currently any licensing and usage of this software is governed under the regulations of the afore mentioned parties

#Author:	Keaghan Brown (3687524) - MSc Bioinformatics Candidate (3687524@myuwc.ac.za)
#Author:	Ruben Cloete (Supervisor) - Lecturer at South African National Bioinformatics Institute (ruben@sanbi.ac.za)

# please ensure that the AMIA folder has been added to the python PYTHONPATH
import sys, getopt
import os
from amia import FoldXAnalysis
from amia import Mutations
from amia import MutantAnalysis
import logging
import __main__
q = MutantAnalysis()
p = Mutations()
m = FoldXAnalysis()

def main (argv):
	inputfile = ''
	try:
		opts, args = getopt.getopt(argv, "hi:",['ifile='])
	except getopt.GetoptError:
		print('amia2.py -i <inputfile>')
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt == '-h':
			print('test.py -i <inputfile> -o <outputfile>')
			sys.exit()
	param_data = open(inputfile).read().split('\n')
	for i in param_data:
		if 'mutation_file' in i:
			mutantion_list_path = i.split('=')[1]
		if 'model_file' in i:
			model_file_path =  i.split('=')[1]
		if 'mutant_storage' in i:
			structure_storage_path = i.split('=')[1]
		if 'introduction_method' in i:
			introduction_mode = i.split('=')[1]
		if 'contacts_output' in i:
			contacts_ouput_folder = i.split('=')[1]
		if 'contact_distance' in i:
			contact_distance  = i.split('=')[1]
	logging.basicConfig(format = '%(asctime)s - %(levelname)s - %(message)s', filename = (os.getcwd() + '/AMIA.log'), filemode='w',  level=logging.DEBUG)
	return mutantion_list_path, model_file_path, structure_storage_path, introduction_mode, contacts_ouput_folder, contact_distance 

			
if __name__ == "__main__":
   	main(sys.argv[1:])

def __init__():
	three_letter ={'V': 'VAL', 'I': 'ILE', 'L': 'LEU', 'E': 'GLU', 'Q': 'GLN', 'D': 'ASP', 'N': 'ASN', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE', 'Y': 'TYR', 'R': 'ARG', 'K': 'LYS', 'S': 'SER', 'T': 'THR', 'M': 'MET', 'A': 'ALA', 'G': 'GLY', 'P': 'PRO', 'C': 'CYS'}
	pymol.cmd.reinitialize()


processing_param_file = main(sys.argv[1:])
p.process_mutant_list(processing_param_file[0])
modes = processing_param_file[3]

"""
if "single" == modes:
    #p.single_mutations(processing_param_file[1], processing_param_file[2])
    #q.single_residue_contacts(processing_param_file[2],processing_param_file[1], processing_param_file[5],processing_param_file[4])
    #q.single_ligand_contacts_calc(processing_param_file[2], processing_param_file[1], processing_param_file[4])
    #m.foldx_stability(processing_param_file[2], processing_param_file[1], processing_param_file[4])
    """
if "multiple" == modes:
    p.multi_mutations(processing_param_file[1], processing_param_file[2])
    p.foldx_opti(processing_param_file[2])
    q.multiple_residue_contacts(processing_param_file[2],processing_param_file[1], processing_param_file[5], processing_param_file[4])
    #q.multi_ligand_contacts_calc(processing_param_file[2], processing_param_file[1], processing_param_file[4])
    #m.foldx_stability(processing_param_file[2], processing_param_file[1], processing_param_file[4])
