# Project Title: "Automated computational workflow to prioritize potential resistance variants identified in HIV Integrase Subtype C and CRF02_AG" 
# This script is developed for the fufuillment for Masters at the South African National Bioinformatics Institute at the University of the Western Cape
# The project is funded by the Poliomyelitis Research Foundation and the UWC Ada & Bertie Levenstein Bursary Programme
# Currently any licensing and usage of this software is governed under the regulations of the afore mentioned parties

#Author:	Keaghan Brown (3687524) - MSc Bioinformatics Candidate (3687524@myuwc.ac.za)
#Author:	Ruben Cloete (Supervisor) - Lecturer at South African National Bioinformatics Institute (ruben@sanbi.ac.za)

import sys, getopt
import os
from tabulate import tabulate
import logging
import __main__

class FoldXAnalysis:
    
    def foldx_stability(self, pdb_path, model_path, output_path, foldx_exe):
        cwd = os.getcwd()
        pdb_path = pdb_path
        model_pathway = os.path.dirname(model_path)
        model = os.path.basename(model_path)
        output_path = output_path
        model_foldx_stability = str(foldx_exe) + ' --command=Stability  --pdb-dir=' + model_pathway + ' --pdb=' + model + ' --output-dir=' + str(output_path) 
        os.chdir(model_pathway)
        os.system(model_foldx_stability)
        energy_table_dict = {"Mutant/ Variant PDB" : [] , "Model ΔG" : [], "Mutant ΔG" : [], "ΔΔG" : []}
        mutant_structures = []
        for file in os.listdir(pdb_path):
            os.chdir(pdb_path)
            if file.endswith('.pdb'):
                mutant_structures.append(file)
                energy_table_dict.setdefault("Mutant/ Variant PDB", []).append(file)
                foldx_stability = str(foldx_exe) + ' --command=Stability  --pdb=' + file  + ' --output-dir=' + output_path
                os.system(foldx_stability)
        os.chdir(output_path)
        model_energy = ''
        for file in os.listdir(output_path):
            if model.split('.')[0] in file:
                model_file = open(file)
                model_data = model_file.read()
                model_energy = model_data.split('\t')[1]
        for file in os.listdir(output_path):
            for mutation in mutant_structures:
                if mutation.split('.')[0] in file:
                    mutant_file = open(file)
                    mutant_data = mutant_file.read()
                    mutant_energy = mutant_data.split('\t')[1]
                    energy_table_dict.setdefault('Mutant ΔG', []).append(mutant_energy)
                    energy_table_dict.setdefault('Model ΔG', []).append(model_energy)
                    ΔΔG =  float(model_energy) - float(mutant_energy)
                    energy_table_dict.setdefault('ΔΔG', []).append(round(ΔΔG, 3))
                    energy_file = open(output_path + '/ΔΔG_output.txt', 'w')
                    energy_file.write(tabulate(energy_table_dict, headers='keys'))
                    energy_file.close()
        
        
