# Project Title: "Automated computational workflow to prioritize potential resistance variants identified in HIV Integrase Subtype C and CRF02_AG" 
# This script is developed for the fufuillment for Masters at the South African National Bioinformatics Institute at the University of the Western Cape
# The project is funded by the Poliomyelitis Research Foundation and the UWC Ada & Bertie Levenstein Bursary Programme
# Currently any licensing and usage of this software is governed under the regulations of the afore mentioned parties

#Author:	Keaghan Brown (3687524) - MSc Bioinformatics Candidate (3687524@myuwc.ac.za)
#Author:	Ruben Cloete (Supervisor) - Lecturer at South African National Bioinformatics Institute (ruben@sanbi.ac.za)

# This script completes the Second Objective for The Automated Computational Workflow to Prioritize Potential Resistance Variants Identified in HIV Integrase Subtype C and CRF02_AG

import pymol
import __main__
__main__.pymol_argv = [ 'pymol', '-Qc']
import logging
import os
from tabulate import tabulate
tabulate.PRESERVE_WHITESPACE = True
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import ntpath
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
import re
from get_raw_distances import *
import itertools

class MutantAnalysis:
	
    def __init__(self):
        self.three_letter ={'V': 'VAL', 'I': 'ILE', 'L': 'LEU', 'E': 'GLU', 'Q': 'GLN', 'D': 'ASP', 'N': 'ASN', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE', 'Y': 'TYR', 'R': 'ARG', 'K': 'LYS', 'S': 'SER', 'T': 'THR', 'M': 'MET', 'A': 'ALA', 'G': 'GLY', 'P': 'PRO', 'C': 'CYS'}
        pymol.cmd.reinitialize('everything')
        self.three_to_one ={'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    def pre_processing(self, file):
        self.file = file
        if self.file.endswith('.pdb'):
            target_mutations = self.file.split(':')[1]
            parser = PDBParser(PERMISSIVE=1)
            structure_id = ntpath.basename(self.file).split('.')[0]
            filename = ntpath.basename(self.file)
            structure = parser.get_structure(structure_id, filename)
            chain_id = []
            chains = structure[0]
            for chain in chains:
                chain_id.append(chain.get_id())
            ppb=PPBuilder()
            ppb.build_peptides(structure)
            sequences = []
            for pp in ppb.build_peptides(structure):
                sequence= pp.get_sequence()
                sequences.append(sequence)
            largest_seq =  ''
            for i in sequences:
                if len(i) > len(largest_seq):
                    largest_seq = i
            mut_resi_pos =  structure_id.split(':')[1]
            if '_' in mut_resi_pos:
                mut_resi_pos = mut_resi_pos.split('_')
            return structure_id, largest_seq, mut_resi_pos,  chain_id, structure
    # This function processes the Mutant files and return a tuple contauning the structure name, the longest sequence present within the PDB file, the mutations present within the structure and the biopython structure objectdef pre_processing(self, file):
    
    def pre_processing2(self, file):
        self.file = file
        if self.file.endswith('.pdb'):
            parser = PDBParser(PERMISSIVE=1)
            structure_id = ntpath.basename(self.file).split('.')[0]
            filename = ntpath.basename(self.file)
            structure = parser.get_structure(structure_id, filename)
            chain_id = []
            chains = structure[0]
            for chain in chains:
                chain_id.append(chain.get_id())
            ppb=PPBuilder()
            ppb.build_peptides(structure)
            sequences = []
            for pp in ppb.build_peptides(structure):
                sequence= pp.get_sequence()
                sequences.append(sequence)
            largest_seq =  ''
            for i in sequences:
                if len(i) > len(largest_seq):
                    largest_seq = i
            return structure_id, largest_seq,  chain_id, structure
    # This function processes the Mutant files and return a tuple contauning the structure name, the longest sequence present within the PDB file, the mutations present within the structure and the biopython structure object

    def contact_calculator2(self, mutant, target_residue,  resi_pos, chains, cutoff):
        self.mutant2 = mutant
        self.chains = chains
        self.resi_pos = resi_pos
        self.target_residue = target_residue
        self.cutoff = cutoff
        contacts = {}
        for i in self.chains:
            surrounding_residues = ''
            pymol.cmd.select(self.target_residue + '_' + str(i),  'resn ' + str(self.three_letter[self.target_residue]) + ' and resi ' + str(self.resi_pos) + ' and chain ' + str(i) + ' around 6 and not resn A+C+G+T')
            seq = pymol.cmd.get_fastastr(str(self.target_residue + '_' + str(i)))
            seq3 = seq.split('\n')
            if len(seq3)>1:
                surrounding_residues = seq3[1]
            surrounding_residues2 = []
            for j in surrounding_residues:
                if j not in surrounding_residues2:
                    surrounding_residues2.append(j)
            if len(surrounding_residues2) > 0:
                for resi in surrounding_residues2:
                    pymol.cmd.select(str(resi) +  '_' + str(i),  'resn ' + str(self.three_letter[resi]) + ' and chain ' + str(i))
                    pymol.cmd.distance('dist_'+ self.mutant2 + str(i) + '_' + str(resi), str(resi) + '_' + str(i), self.mutant2, '3.6', mode='2')
                    D = get_raw_distances('dist_'+ self.mutant2 + str(i) + '_' +  str(resi))
                    contacts[str(i) + '|' + str(resi)] = (len(D))
                    pymol.cmd.delete(str(resi) +  '_' + str(i))
                    pymol.cmd.delete('dist_'+ self.mutant2 + str(i) + '_' + str(resi))
        return contacts
        
    def pdb_ligands(self, file):
        self.file = file
        ligand_res = []
        DNA = ['A', 'C', 'G', 'T']
        structure = ''
        if ':' in self.file:
            ligands = self.pre_processing(self.file)
            structure = ligands[4]
        elif ':' not in self.file:
            ligands = self.pre_processing2(self.file)
            structure = ligands[3]
        for model in structure:
            residues = model.get_residues()
            for residue in residues:
                residue = residue.get_resname()
                if residue not in self.three_to_one:    
                    if residue not in DNA and not residue in ligand_res:
                        ligand_res.append(residue)
        return(ligand_res, DNA)
                    
    def ligand_contacts(self, file):
        self.file = file 
        contacts = {}
        non_prot_data = self.pdb_ligands(self.file)
        DNA = non_prot_data[1]
        DNA_seq =''
        for nucleotide in  DNA:
            DNA_seq += nucleotide + '+'
        if DNA_seq.endswith('+'):
            DNA_seq = DNA_seq[0:-1]
        ligands = non_prot_data[0]
        pymol.cmd.select('DNA_' +  str(DNA_seq), 'resn ' + str(DNA_seq))
        pymol.cmd.select('Around_DNA_' +  str(DNA_seq), 'resn ' + str(DNA_seq) + ' around 6')
        pymol.cmd.distance('Contacts', 'Around_DNA_' +  str(DNA_seq) , 'DNA_' +  str(DNA_seq), '3.6', mode = '2')
        D3 = get_raw_distances('Contacts',)
        contacts['contact_@_DNA'] = (len(D3))
        for lig in ligands:
            if ':' in file:
                for chain in self.pre_processing(self.file)[3]:
                    pymol.cmd.select('ligand_' + str(lig) + str(chain), 'resn ' + str(lig) + ' and chain ' + str(chain))
                    pymol.cmd.select('ligand2_' + str(lig) + str(chain), 'resn ' + str(lig) +  ' and chain ' + str(chain)+  '  around 3.5')
                    contact_length = pymol.cmd.distance('contacts_' + 'ligand_' + str(lig) + str(chain) + 'ligand2_' + str(lig) + str(chain), 'ligand_' + str(lig) + str(chain), 'ligand2_' + str(lig) + str(chain), '3.5', mode=2)
                    if contact_length > 0:
                        D4 = get_raw_distances('contacts_' + 'ligand_' + str(lig) + str(chain) + 'ligand2_' + str(lig) + str(chain))
                        contacts['contact_@_' + str(lig) + '_' +  str(chain)] = (len(D4))
            if ':' not in file:
                for chain in self.pre_processing2(self.file)[2]:
                    pymol.cmd.select('ligand_' + str(lig) + str(chain), 'resn ' + str(lig) + ' and chain ' + str(chain))
                    pymol.cmd.select('ligand2_' + str(lig) + str(chain), 'resn ' + str(lig) +  ' and chain ' + str(chain)+  '  around 3.5')
                    contact_length = pymol.cmd.distance('contacts_' + 'ligand_' + str(lig) + str(chain) + 'ligand2_' + str(lig) + str(chain), 'ligand_' + str(lig) + str(chain), 'ligand2_' + str(lig) + str(chain), '3.5', mode=2)
                    if contact_length > 0:
                        D4 = get_raw_distances('contacts_' + 'ligand_' + str(lig) + str(chain) + 'ligand2_' + str(lig) + str(chain))
                        contacts['contact_@_' + str(lig) + '_' +  str(chain)] = (len(D4))
        return contacts
        

    def single_residue_contacts(self, mutant_path, model, cut_off, output):
        pymol.finish_launching()
        self.mutant_path = mutant_path
        self.output = output
        self.model = model
        self.cut_off = cut_off
        table_dict2 = {"Mutant" : [], "Wild - Type Contacts" : [],  "Mutant Contacts" : [],  "No. of Wild - Type Contacts" : [], "No. of Mutant Contacts" : []}
        data_file = open(output + '/single_mutant_contacts.txt', 'w')
        data_file.write(tabulate(table_dict2, headers="keys", numalign = "center", stralign = "center"))
        logging.basicConfig(format = '%(asctime)s - %(levelname)s - %(message)s', filename = (output+ '/single_resi_conts.log'), filemode='w',  level=logging.DEBUG)
        logging.info("RUNNING SINGLE MUTANT CONTACT ANALYSIS....")
        for file in os.listdir(self.mutant_path):
            os.chdir(self.mutant_path)
            if file.endswith('.pdb'):
                logging.info("Structure " + file + " has been loaded for intra-residue contact analysis")
                analysis = self.pre_processing(file)
                resi_pos = analysis[2][1:-1]
                mutated_residue = analysis[2][-1]
                initial_residue = analysis[2][0]
                table_dict = {"Mutant" : [], "Wild - Type Contacts" : [],  "Mutant Contacts" : [],  "No. of Wild - Type Contacts"  : [], "No. of Mutant Contacts" : []}
                table_dict.setdefault("Mutant", []).append(analysis[2])
                pymol.cmd.reinitialize('everything')
                pymol.cmd.load(file)
                pymol.cmd.select('mutant_pos_'  + str(analysis[2][1:-1]) , ' resn ' + str(self.three_letter[analysis[2][-1]]) + ' and resi ' +  str(analysis[2][1:-1]))
                mutant_contacts = self.contact_calculator2('mutant_pos_'  + str(analysis[2][1:-1]), analysis[2][-1], analysis[2][1:-1], analysis[3], '5')
                for contact in mutant_contacts:
                    if mutant_contacts[contact] > 0:
                        table_dict.setdefault("Mutant Contacts", []).append(contact)
                        table_dict.setdefault('No. of Mutant Contacts', []).append(str(mutant_contacts[contact]))
                        os.chdir(os.path.dirname(self.model))
                pymol.cmd.reinitialize('everything')
                pymol.cmd.load(self.model)
                pymol.cmd.select('model_pos_'  + str(analysis[2][1:-1]) , ' resn ' + str(self.three_letter[analysis[2][0]]) + ' and resi ' +  str(analysis[2][1:-1]))
                model_contacts = self.contact_calculator2('model_pos_'  + str(analysis[2][1:-1]), analysis[2][0],  analysis[2][1:-1], analysis[3], '5')
                for contact2 in model_contacts:
                    if model_contacts[contact2] > 0:
                        table_dict.setdefault('Wild - Type Contacts', []).append(contact2)
                        table_dict.setdefault('No. of Wild - Type Contacts', []).append(str(model_contacts[contact2]))
                for i in table_dict["Wild - Type Contacts"]:
                    if i not in table_dict["Mutant Contacts"] and i != '0':
                        pos = table_dict["Wild - Type Contacts"].index(i)
                        table_dict["Mutant Contacts"].insert(pos, '0')
                        table_dict["No. of Mutant Contacts"].insert(pos, '0')
                for j in table_dict["Mutant Contacts"]:
                    if j != '0' and j not in table_dict["Wild - Type Contacts"]  :
                        pos = table_dict["Mutant Contacts"].index(j)
                        table_dict["Wild - Type Contacts"].insert(pos, '0')
                        table_dict["No. of Wild - Type Contacts"].insert(pos, '0')
                new_headers = []
                for key in table_dict:
                    new_key = ''
                    for i in key:
                        new_key += ' '
                    new_headers.append(new_key)
                data_file = open(output + '/single_mutant_contacts.txt', 'a')
                data_file.write(tabulate(table_dict, headers= new_headers, numalign = "center", stralign = "center", tablefmt = "plain")+ '\n' )
                data_file.close()

    def multiple_residue_contacts(self, mutant_folder, model, cut_off, output):
        pymol.finish_launching()
        pymol.cmd.reinitialize('everything')
        logging.basicConfig(format = '%(asctime)s - %(levelname)s - %(message)s', filename = (output + '/multi_resi_conts.log'), filemode='w',  level=logging.DEBUG)
        logging.info("RUNNING MULTIPLE MUTANT CONTACT ANALYSIS....")
        self.mutant_folder_path = mutant_folder
        self.model = model
        self.cut_off = cut_off
        model_dir = os.path.dirname(self.model)
        mutant_dir = self.mutant_folder_path
        table_dict2 = {"Mutant" : [], "Wild - Type Contacts" : [],  "Mutant Contacts" : [],  "No. of Wild - Type Contacts" : [], "No. of Mutant Contacts" : []}
        data_file2 = open(output + '/multiple_mutant_contacts.txt', 'w')
        data_file2.write(tabulate(table_dict2, headers="keys", numalign = "center", stralign = "center") + '\n')
        for file in os.listdir(mutant_dir):
            os.chdir(mutant_dir)
            if file.endswith('.pdb'):
                logging.info("Structure " + file + " has been loaded for intra-residue contact analysis")
                analysis = self.pre_processing(file)
                if type(analysis[2]) is str:
                    table_dict = {"Mutant" : [], "Wild - Type Contacts" : [],  "Mutant Contacts" : [],  "No. of Wild - Type Contacts"  : [], "No. of Mutant Contacts" : []}
                    resi_pos = analysis[2][1:-1]
                    mutated_residue = analysis[2][-1]
                    initial_residue = analysis[2][0]
                    table_dict.setdefault("Mutant", []).append(analysis[2])
                    pymol.cmd.reinitialize('everything')
                    pymol.cmd.load(file)
                    pymol.cmd.select('mutant_pos_'  + str(analysis[2][1:-1]) , ' resn ' + str(self.three_letter[analysis[2][-1]]) + ' and resi ' +  str(analysis[2][1:-1]))
                    mutant_contacts = self.contact_calculator2('mutant_pos_'  + str(analysis[2][1:-1]), analysis[2][-1], analysis[2][1:-1], analysis[3], '5')
                    for contact in mutant_contacts:
                        if mutant_contacts[contact] > 0:
                            table_dict.setdefault("Mutant Contacts", []).append(contact)
                            table_dict.setdefault('No. of Mutant Contacts', []).append(str(mutant_contacts[contact]))
                            os.chdir(os.path.dirname(self.model))
                    pymol.cmd.reinitialize('everything')
                    pymol.cmd.load(self.model)
                    pymol.cmd.select('model_pos_'  + str(analysis[2][1:-1]) , ' resn ' + str(self.three_letter[analysis[2][0]]) + ' and resi ' +  str(analysis[2][1:-1]))
                    model_contacts = self.contact_calculator2('model_pos_'  + str(analysis[2][1:-1]), analysis[2][0],  analysis[2][1:-1], analysis[3], '5')
                    for contact2 in model_contacts:
                        if model_contacts[contact2] > 0:
                            table_dict.setdefault('Wild - Type Contacts', []).append(contact2)
                            table_dict.setdefault('No. of Wild - Type Contacts', []).append(str(model_contacts[contact2]))
                    for i in table_dict["Wild - Type Contacts"]:
                        if i not in table_dict["Mutant Contacts"] and i != '0':
                            pos = table_dict["Wild - Type Contacts"].index(i)
                            table_dict["Mutant Contacts"].insert(pos, '0')
                            table_dict["No. of Mutant Contacts"].insert(pos, '0')
                    for j in table_dict["Mutant Contacts"]:
                        if j != '0' and j not in table_dict["Wild - Type Contacts"]  :
                            pos = table_dict["Mutant Contacts"].index(j)
                            table_dict["Wild - Type Contacts"].insert(pos, '0')
                            table_dict["No. of Wild - Type Contacts"].insert(pos, '0')
                    new_headers = []
                    for key in table_dict:
                        new_key = ''
                        for i in key:
                            new_key += ' '
                        new_headers.append(new_key)
                    data_file2 = open(output + '/multiple_mutant_contacts.txt',  'a')
                    data_file2.write(tabulate(table_dict, headers= new_headers, numalign = "center", stralign = "center", tablefmt = "plain")+ '\n' )
                elif type(analysis[2]) is list:
                    for mutation in analysis[2]:
                        os.chdir(mutant_dir)
                        table_dict = {"Mutant" : [], "Wild - Type Contacts" : [],  "Mutant Contacts" : [],  "No. of Wild - Type Contacts"  : [], "No. of Mutant Contacts" : []}
                        resi_pos = mutation[1:-1]
                        mutated_residue = mutation[-1]
                        initial_residue = mutation[0]
                        table_dict.setdefault("Mutant", []).append(mutation)
                        pymol.cmd.reinitialize('everything')
                        pymol.cmd.load(file)
                        pymol.cmd.select('mutant_pos_'  + str(resi_pos) , ' resn ' + str(self.three_letter[mutated_residue]) + ' and resi ' +  str(resi_pos))
                        mutant_contacts = self.contact_calculator2('mutant_pos_'  + str(resi_pos),mutated_residue, resi_pos, analysis[3], '5')
                        for contact in mutant_contacts:
                            if mutant_contacts[contact] > 0:
                                table_dict.setdefault("Mutant Contacts", []).append(contact)
                                table_dict.setdefault('No. of Mutant Contacts', []).append(str(mutant_contacts[contact]))
                                os.chdir(os.path.dirname(self.model))
                        pymol.cmd.reinitialize('everything')
                        pymol.cmd.load(self.model)
                        pymol.cmd.select('model_pos_'  + str(mutation[1:-1]) , ' resn ' + str(self.three_letter[mutation[0]]) + ' and resi ' +  str(mutation[1:-1]))
                        model_contacts = self.contact_calculator2('model_pos_'  + str(mutation[1:-1]), mutation[0],  mutation[1:-1], analysis[3], '5')
                        new_mut = []
                        new_wt = []
                        for contact2 in model_contacts:
                            if model_contacts[contact2] > 0:
                                table_dict.setdefault('Wild - Type Contacts', []).append(contact2)
                                table_dict.setdefault('No. of Wild - Type Contacts', []).append(str(model_contacts[contact2]))
                        for i, j in itertools.zip_longest(table_dict["Wild - Type Contacts"], table_dict["Mutant Contacts"]):
                            if i in table_dict["Mutant Contacts"] and  i != None:
                                new_mut.append(i)
                                new_wt.append(i)
                            elif i not in table_dict ["Mutant Contacts"] and i != None:
                                new_mut.append('0')
                                new_wt.append(i)
                            if j not in table_dict["Wild - Type Contacts"] and j != None :
                                new_mut.append(j)
                                new_wt.append('0')
                            elif j not in table_dict["Wild - Type Contacts"] and  j != None:
                                new_wt.append('0')
                                new_mut.append(j)
                        table_dict["Mutant Contacts"] = new_mut
                        table_dict["Wild - Type Contacts"] = new_wt
                        for k in range(len(new_wt)):
                            if new_wt[k] == '0':
                                table_dict["No. of Wild - Type Contacts"].insert(k, '0')
                        for l in range(len(new_mut)):
                            if new_mut[l] == '0':
                                table_dict["No. of Mutant Contacts"].insert(l, '0')
                        new_headers = []
                        for key in table_dict:
                            new_key = ''
                            for i in key:
                                new_key += ' '
                            new_headers.append(new_key)
                        data_file2 = open(output + '/multiple_mutant_contacts.txt',  'a')
                        data_file2.write(tabulate(table_dict, headers= new_headers, numalign = "center", stralign = "center", tablefmt = "plain")+ '\n' )
            data_file2.close() 
                
                
    def multi_ligand_contacts_calc(self, mutant_folder, model, output):
        pymol.cmd.reinitialize('everything')
        self.mutant_folder = mutant_folder
        self.model = model
        self.output = output
        table_dict3 = {"Mutant" : [], "Wild - Type Contacts" : [],  "Mutant Contacts" : [],  "No. of Wild - Type Contacts" : [], "No. of Mutant Contacts" : []}
        data_file3 = open(output + '/multi_mutant_ligand_contacts.txt', 'w')
        data_file3.write(tabulate(table_dict3, headers="keys", numalign = "center", stralign = "center"))
        for file in os.listdir(self.mutant_folder):
            pymol.cmd.reinitialize('everything')
            os.chdir(self.mutant_folder)
            if file.endswith('.pdb'):
                logging.info("Structure " + file + " has been loaded for ligandcontact analysis")
                analysis = self.pre_processing(file)
                table_dict4 = {"Mutant" : [], "Wild - Type Contacts" : [],  "Mutant Contacts" : [],  "No. of Wild - Type Contacts"  : [], "No. of Mutant Contacts" : []}
                if type(analysis[2]) is list:
                    for mutant in analysis[2]:
                        table_dict4.setdefault("Mutant", []).append(mutant)
                elif type(analysis[2]) is str:
                    table_dict4.setdefault("Mutant", []).append(analysis[2])
                pymol.cmd.reinitialize('everything')
                pymol.cmd.load(file)
                mutant_non_residue_contacts = self.ligand_contacts(file)
                for contact3 in mutant_non_residue_contacts :
                    table_dict4.setdefault("Mutant Contacts", []).append(contact3)
                    table_dict4.setdefault('No. of Mutant Contacts', []).append(str(mutant_non_residue_contacts [contact3]))
                os.chdir(os.path.dirname(self.model))
                pymol.cmd.reinitialize('everything')
                pymol.cmd.load(self.model)
                wt_non_residue_contacts = self.ligand_contacts(self.model)
                for contact4 in wt_non_residue_contacts:
                    table_dict4.setdefault("Wild - Type Contacts", []).append(contact4)
                    table_dict4.setdefault('No. of Wild - Type Contacts', []).append(str(wt_non_residue_contacts[contact4]))
                new_headers = []
                for key in table_dict3:
                    new_key = ''
                    for i in key:
                        new_key += ' '
                    new_headers.append(new_key)
                data_file3 = open(output + '/multi_mutant_ligand_contacts.txt', 'a')
                data_file3.write(tabulate(table_dict4, headers= new_headers, numalign = "center", stralign = "center", tablefmt = "plain")+ '\n' )
            data_file3.close()
                
    def single_ligand_contacts_calc(self, mutant_folder, model, output):
        self.mutant_folder = mutant_folder
        self.model = model
        self.output = output
        table_dict3 = {"Mutant" : [], "Wild - Type Contacts" : [],  "Mutant Contacts" : [],  "No. of Wild - Type Contacts" : [], "No. of Mutant Contacts" : []}
        data_file3 = open(output + '/single_mutant_ligand_contacts.txt', 'w')
        data_file3.write(tabulate(table_dict3, headers="keys", numalign = "center", stralign = "center"))
        for file in os.listdir(self.mutant_folder):
            logging.info("Structure " + file + " has been loaded for ligand contact analysis")
            os.chdir(self.mutant_folder)
            if file.endswith('.pdb'):
                analysis = self.pre_processing(file)
                table_dict4 = {"Mutant" : [], "Wild - Type Contacts" : [],  "Mutant Contacts" : [],  "No. of Wild - Type Contacts"  : [], "No. of Mutant Contacts" : []}
                table_dict4.setdefault("Mutant", []).append(analysis[2])
                pymol.cmd.reinitialize('everything')
                pymol.cmd.load(file)
                mutant_non_residue_contacts = self.ligand_contacts(file)
                for contact3 in mutant_non_residue_contacts :
                    table_dict4.setdefault("Mutant Contacts", []).append(contact3)
                    table_dict4.setdefault('No. of Mutant Contacts', []).append(str(mutant_non_residue_contacts [contact3]))
                os.chdir(os.path.dirname(self.model))
                pymol.cmd.reinitialize('everything')
                pymol.cmd.load(self.model)
                wt_non_residue_contacts = self.ligand_contacts(self.model)
                for contact4 in wt_non_residue_contacts:
                    table_dict4.setdefault("Wild - Type Contacts", []).append(contact4)
                    table_dict4.setdefault('No. of Wild - Type Contacts', []).append(str(wt_non_residue_contacts[contact4]))
                pymol.cmd.reinitialize()
                pymol.cmd.load(self.model)
                new_headers = []
                for key in table_dict3:
                    new_key = ''
                    for i in key:
                        new_key += ' '
                    new_headers.append(new_key)
                data_file3 = open(output + '/single_mutant_ligand_contacts.txt', 'a')
                data_file3.write(tabulate(table_dict4, headers= new_headers, numalign = "center", stralign = "center", tablefmt = "plain")+ '\n' )
                data_file3.close()