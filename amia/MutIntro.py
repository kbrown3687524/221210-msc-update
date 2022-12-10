# Project Title: "Automated computational workflow to prioritize potential resistance variants identified in HIV Integrase Subtype C and CRF02_AG" 
# This script is developed for the fufuillment for Masters at the South African National Bioinformatics Institute at the University of the Western Cape
# The project is funded by the Poliomyelitis Research Foundation and the UWC Ada & Bertie Levenstein Bursary Programme
# Currently any licensing and usage of this software is governed under the regulations of the afore mentioned parties

#Author:	Keaghan Brown (3687524) - MSc Bioinformatics Candidate (3687524@myuwc.ac.za)
#Author:	Ruben Cloete (Supervisor) - Lecturer at South African National Bioinformatics Institute (ruben@sanbi.ac.za)

# This script completes the first objective for The Automated Computational Workflow to Prioritize Potential Resistance Variants Identified in HIV Integrase Subtype C and CRF02_AG

import pymol
from Bio.PDB.PDBParser import PDBParser 
from Bio.PDB.Polypeptide import PPBuilder
import ntpath
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
import re
import logging
import os
from pathlib import Path
import __main__
__main__.pymol_argv = [ 'pymol', '-Qc' ]
import pymol
pymol.finish_launching()
from pymol import cmd

class Mutations:

    def __init__(self):
        self.single_list = {}
        self.multi_list = {}
        self.three_letter ={'V': 'VAL', 'I': 'ILE', 'L': 'LEU', 'E': 'GLU', 'Q': 'GLN', 'D': 'ASP', 'N': 'ASN', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE', 'Y': 'TYR', 'R': 'ARG', 'K': 'LYS', 'S': 'SER', 'T': 'THR', 'M': 'MET', 'A': 'ALA', 'G': 'GLY', 'P': 'PRO', 'C': 'CYS'}
        
    def process_mutant_list(self, mutant_list):
        mut_path = os.path.dirname(mutant_list)
        os.chdir(mut_path)
        if not re.search(r'.txt', mutant_list):
            logging.error("A suitable file formt not supplied for mutations")
            raise Exception("Mutation file not supplied, please view log file")
        if os.path.basename(mutant_list) in os.listdir(mut_path):
            data = open(mutant_list)
            file_data = data.read()
            patient_data = file_data.split('$')
            patient_data.remove('') 
            for patient in patient_data:
                data = patient.split('^^^')
                source = data[0].strip()
                mutants = data[1].split('\n')
                while '' in mutants:
                    mutants.remove('')
                self.single_list[source] = mutants
                self.multi_list[source] = mutants
        else:
            logging.error("Could not process mutation list due to incorrect format, please ensure format matches example file")
        return


	# The script above, attempts to retrieve a user specified list of mutations in a standard .txt file. The format of the file is as follows: '$' denotes the source of the mutations and for which drug they correlate to, '^^^'  denotes the start of the mutation list. Once each subset has been seperated, they are added to a dictionary which will be used by other functions for mutation introduction
    def single_mutations(self, model, output):
        model_path = os.path.dirname(model)
        os.chdir(model_path)
        logging.basicConfig(format = '%(asctime)s - %(levelname)s - %(message)s', filename = (output+ '/single_mutant_intro.log'), filemode='w',  level=logging.DEBUG)
        logging.info("SINGLE MUTATION INTRODUCTION METHOD IS RUNNING....")
        if not re.search(r'.pdb',model):
            logging.error("A suitable PDB file was not supplied for mutation introduction")
            raise Exception("PDB file not supplied, please view log file")
        if os.path.basename(model) in os.listdir(model_path):
            structure_id = os.path.basename(model)
            file_name = os.path.basename(model)
            parser = PDBParser(PERMISSIVE = 1)
            structure = parser.get_structure(structure_id, file_name)
            ppbuilder = PPBuilder()
            sequences = []
            for pp in ppbuilder.build_peptides(structure):
                sequence= pp.get_sequence()
                sequences.append(sequence)
            chain_start_pos = pp[0].get_id()[1]
            largest_seq =  ''
            for i in  sequences:
                if len(i) > len(largest_seq):
                    largest_seq = i
            for key in self.single_list:
                for mutation in self.single_list[key]:
                    print('Attempting to introduce ' + mutation + '  mutation....')
                    initial_residue = mutation[0]
                    mutated_residue = mutation[len(mutation)-1]
                    residue_pos = mutation[1:len(mutation)-1]
                    if int(residue_pos) <= int(len(largest_seq)):
                        cmd.reinitialize()
                        cmd.load(model)
                        for i in cmd.get_chains(ntpath.basename(model).split('.')[0]):
                            cmd.select('mutant_' + mutation, 'resn ' + str(self.three_letter[initial_residue]) + ' and resi ' + residue_pos)
                            cmd.wizard('mutagenesis')
                            cmd.refresh_wizard()
                            cmd.get_wizard().do_select("/" + 'mutant_' + mutation + "//" + i  + "/" + residue_pos + "/")
                            cmd.get_wizard().set_mode(self.three_letter[mutated_residue])
                            cmd.get_wizard().apply()
                            print(mutation + '  mutation successfully introduced')
                            print('Attempting to energy minimize ' + str(mutation) + '.pdb structure....')
                        cmd.select('mutation',  'resn ' + str(self.three_letter[initial_residue]) + ' and resi ' + str(residue_pos) + ' around 3')
                        cmd.set_wizard()
                        os.chdir(output)
                        logging.info("Mutation " + mutation + " PDB was successfully generated from the WT structure!")
                        cmd.save(key + str(mutation) + '.pdb')
                        print(str(mutation) + '.pdb structure saved to ' + output)
                    elif int(residue_pos) > int(len(largest_seq)):
                        logging.error("Mutation " + mutation + " position exceeds chain length")
        cmd.set_wizard("done")
        cmd.reinitialize()  
        return
                    
# The script above, attempts to introduce the previosly obtained mutations, individually into the supplied PDB structure by reloading the WT structure with each mutant. Each mutation is processed into the initial residue, residue poition and mutated residue. The script runs through each chain in the PDB and specifically selects the initial residue and position and stores it as an object, which ensures that only residues are mutated and not DNA, ligands or ions. The script then calls on the mutagenesis wizard to mutate the selected residues in the object from the initial residue to the mutant. As PyMOL automatically selects the rotamer with the least sterich clashes, we can employ the programs built-in FF to minimize the rersidue and surrounding region (3A) afterwhich the script generates a mutated PDB in a user specified folder

    def multi_mutations(self, model, output):
        model_path = os.path.dirname(model)
        os.chdir(model_path)
        multi_mutant_list = ''
        cmd.reinitialize()
        logging.basicConfig(format = '%(asctime)s - %(levelname)s - %(message)s', filename = (output+ '/multi_mutant_intro.log'), filemode='w',  level=logging.DEBUG)
        logging.info("MULTIPLE MUTATION INTRODUCTION METHOD IS RUNNING....")
        if not re.search(r'.pdb',model):
            logging.error("A suitable PDB file was not supplied for mutation introduction")
            raise Exception("PDB file not supplied, please view log file")
        if os.path.basename(model) in os.listdir(model_path):
            structure_id = os.path.basename(model)
            file_name = os.path.basename(model)
            parser = PDBParser(PERMISSIVE = 1)
            structure = parser.get_structure(structure_id, file_name)
            ppbuilder = PPBuilder()
            sequences = []
            for pp in ppbuilder.build_peptides(structure):
                sequence= pp.get_sequence()
                sequences.append(sequence)
            chain_start_pos = pp[0].get_id()[1]
            largest_seq =  ''
            for i in  sequences:
                if len(i) > len(largest_seq):
                    largest_seq = i
            for key in self.multi_list:
                mutant_list = ''
                cmd.reinitialize()
                cmd.load(model)
                for mutation in self.multi_list[key]:
                    mutant_list += mutation + '_'
                if mutant_list.endswith('_'):
                    mutant_list = mutant_list[:-1]
                for mutation in self.multi_list[key]:
                    initial_residue = mutation[0]
                    mutated_residue = mutation[len(mutation)-1]
                    residue_pos = mutation[1:len(mutation)-1]
                    if int(residue_pos) <= int(len(largest_seq)):
                        for i in cmd.get_chains(ntpath.basename(model).split('.')[0]):
                            cmd.select('mutant_' + mutation, 'resn ' + str(self.three_letter[initial_residue]) + ' and resi ' + residue_pos)
                            cmd.wizard('mutagenesis')
                            cmd.refresh_wizard()
                            cmd.get_wizard().do_select("/" + 'mutant_' + mutation + "//" + i  + "/" + residue_pos + "/")
                            cmd.get_wizard().set_mode(self.three_letter[mutated_residue])
                            cmd.get_wizard().apply()
                            print(mutation + '  mutation successfully introduced')
                            print('Attempting to energy minimize ' + str(mutation) + '.pdb structure....')
                        cmd.select('mutation',  'resn ' + str(self.three_letter[initial_residue]) + ' and resi ' + str(residue_pos) + ' around 3')
                        cmd.set_wizard()
                        os.chdir(output)
                    elif int(residue_pos) > int(len(largest_seq)):
                        mutant_list = mutant_list.replace(mutation, '')
                        if mutant_list.endswith('_'):
                            mutant_list = mutant_list[:-1]
                        self.multi_list[key].remove(mutation)
                        logging.error("Mutation " + mutation + " position exceeds chain length")
                logging.info("Mutations " + str(self.multi_list[key])+ " PDB was successfully generated from the WT structure!")
                cmd.save(key + mutant_list + '.pdb')
                print(key + mutant_list +  '.pdb structure saved to ' + output)
        cmd.set_wizard("done")
        cmd.reinitialize()
        return
        
    def foldx_opti(self, output):
        for file in os.listdir(output):
            os.chdir(output)
            if file.endswith('.pdb'):
                foldx_optimization = '/home/akuma/Desktop/App/FoldX/foldx --command=Optimize  --pdb=' + file  + ' --output-dir=' + output + '--output-file=' + file
                os.system(foldx_optimization)
            file.replace(':','_')
        return
        
            
# The script above, attempts to introduce the previosly obtained mutations, all at once into the supplied PDB structure by loading the WT structure prior to introducing each mutant. Each mutation is processed into the initial residue, residue poition and mutated residue. The script runs through each chain in the PDB and specifically selects the initial residue and position and stores it as an object, which ensures that only residues are mutated and not DNA, ligands or ions. The script then calls on the mutagenesis wizard to mutate the selected residues in the object from the initial residue to the mutant. As PyMOL automatically selects the rotamer with the least sterich clashes, we can employ the programs built-in FF to minimize the rersidue and surrounding region (3A) afterwhich the script generates a mutated PDB in the user specified folder
