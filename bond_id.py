from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import ntpath
import warnings
import re
from Bio import BiopythonWarning
from center_of_mass import com
from get_raw_distances import *

class BOND_ID:

    def __init__(self):
        cmd.reinitialize()

    def ionic_bonds(self, structure_path):
        three_letter ={'V': 'VAL', 'I': 'ILE', 'L': 'LEU', 'E': 'GLU', 'Q': 'GLN', 'D': 'ASP', 'N': 'ASN', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE', 'Y': 'TYR', 'R': 'ARG', 'K': 'LYS', 'S': 'SER', 'T': 'THR', 'M': 'MET', 'A': 'ALA', 'G': 'GLY', 'P': 'PRO', 'C': 'CYS'}
        self.protein = structure_path
        prot_path = os.path.dirname(os.path.abspath(self.protein))
        os.chdir(prot_path)
        parser = PDBParser(PERMISSIVE=1)
        structure_id = ntpath.basename(self.protein).split('.')[0]
        filename = ntpath.basename(self.protein)
        structure = parser.get_structure(structure_id, filename)
        ppb=PPBuilder()
        ppb.build_peptides(structure)
        cmd.load(self.protein)
        cmd.show('sticks')
        cmd.hide('cartoon')
        sequences = {}
        seq_list = []
        pos_lst = []
        neg_lst=[]
        ionic_bonds = []
        for pp in ppb.build_peptides(structure):
            seq = pp.get_sequence()
            chain_sequence = pp.get_sequence()      
            seq_list.append(chain_sequence)
            chain_start_pos = pp[0].get_id()[1] 
            sequences[seq] = chain_start_pos
        for i in seq_list:
            if sequences[i] != 1:
                counter = sequences[i] 
            else:
                counter = 0
            for j in i:
                counter +=1               
                if j == 'K':
                    for chain in cmd.get_chains(ntpath.basename(self.protein).split('.')[0]):
                        cmd.select('POS_LYS_' + str(counter-1) + '_' + str(chain), 'chain ' + str(chain) + ' and resi ' + str(counter-1) + ' and name NZ')
                        atom_count = cmd.count_atoms('POS_LYS_' + str(counter-1) + '_' + str(chain))
                        if atom_count > 0:
                            com('POS_LYS_' + str(counter-1) + '_' + str(chain))
                            pos_lst.append('POS_LYS_' + str(counter-1) + '_' + str(chain) + '_COM')
        for i in seq_list:
            if sequences[i] != 1:
                counter = sequences[i] 
            else:
                counter = 0
            for j in i:
                counter +=1               
                if j == 'R':
                    for chain in cmd.get_chains(ntpath.basename(self.protein).split('.')[0]):
                        cmd.select('POS_ARG_' + str(counter-1) + '_' + str(chain), 'chain ' + str(chain) + ' and resi ' + str(counter-1) + ' and name NE+NH*')
                        atom_count = cmd.count_atoms('POS_ARG_' + str(counter-1) + '_' + str(chain))
                        if atom_count > 0:
                            com('POS_ARG_' + str(counter-1) + '_' + str(chain))
                            pos_lst.append('POS_ARG_' + str(counter-1) + '_' + str(chain) + '_COM')
        for i in seq_list:
            if sequences[i] != 1:
                counter = sequences[i] 
            else:
                counter = 0
            for j in i:
                counter +=1               
                if j == 'H':
                    for chain in cmd.get_chains(ntpath.basename(self.protein).split('.')[0]):
                        cmd.select('POS_HIS_' + str(counter-1) + '_' + str(chain), 'chain ' + str(chain) + ' and resi ' + str(counter-1) + ' and name NE*+ND*')
                        atom_count = cmd.count_atoms('POS_HIS_' + str(counter-1) + '_' + str(chain))
                        if atom_count > 0:
                            com('POS_HIS_' + str(counter-1) + '_' + str(chain))
                            pos_lst.append('POS_HIS_' + str(counter-1) + '_' + str(chain) + '_COM')
        for i in seq_list:
            if sequences[i] != 1:
                counter = sequences[i] 
            else:
                counter = 0
            for j in i:
                counter +=1               
                if j == 'D':
                    for chain in cmd.get_chains(ntpath.basename(self.protein).split('.')[0]):
                        cmd.select('NEG_ASP_' + str(counter-1) + '_' + str(chain), 'chain ' + str(chain) + ' and resi ' + str(counter-1) + ' and name OD*+OE*')
                        atom_count = cmd.count_atoms('NEG_ASP_' + str(counter-1) + '_' + str(chain))
                        if atom_count > 0:
                            com('NEG_ASP_' + str(counter-1) + '_' + str(chain))
                            neg_lst.append('NEG_ASP_' + str(counter-1) + '_' + str(chain) + '_COM')
        for i in seq_list:
            if sequences[i] != 1:
                counter = sequences[i] 
            else:
                counter = 0
            for j in i:
                counter +=1               
                if j == 'E':
                    for chain in cmd.get_chains(ntpath.basename(self.protein).split('.')[0]):
                        cmd.select('NEG_GLU_' + str(counter-1) + '_' + str(chain), 'chain ' + str(chain) + ' and resi ' + str(counter-1) + ' and name OD*+OE*')
                        atom_count = cmd.count_atoms('NEG_GLU_' + str(counter-1) + '_' + str(chain))
                        if atom_count > 0:
                            com('NEG_GLU_' + str(counter-1) + '_' + str(chain))
                            neg_lst.append('NEG_GLU_' + str(counter-1) + '_' + str(chain) + '_COM')
        for k in pos_lst:
            for l in neg_lst:
                cmd.distance('IONIC_' + str(k)+'_'+ str(l), str(k), str(l), '4', mode=0)
        for k in pos_lst:
            for l in neg_lst:
                cmd.delete(k)
                cmd.delete(l)
            
                    
    def hydrophobic_bonds(self, model):
        hydrophobe_dict = {}
        three_letter ={'V': 'VAL', 'I': 'ILE', 'L': 'LEU', 'E': 'GLU', 'Q': 'GLN', 'D': 'ASP', 'N': 'ASN', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE', 'Y': 'TYR', 'R': 'ARG', 'K': 'LYS', 'S': 'SER', 'T': 'THR', 'M': 'MET', 'A': 'ALA', 'G': 'GLY', 'P': 'PRO', 'C': 'CYS'}
        hydrophobes = ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
        self.protein = model
        prot_path = os.path.dirname(os.path.abspath(self.protein))
        os.chdir(prot_path)
        parser = PDBParser(PERMISSIVE=1)
        structure_id = ntpath.basename(self.protein).split('.')[0]
        filename = ntpath.basename(self.protein)
        structure = parser.get_structure(structure_id, filename)
        ppb=PPBuilder()
        ppb.build_peptides(structure)
        cmd.load(self.protein)
        cmd.show('sticks')
        cmd.hide('cartoon')
        hydrophobe_resi = []
        sequences = {}
        seq_list = []
        pos_lst = []
        neg_lst=[]
        ionic_bonds = []
        for pp in ppb.build_peptides(structure):
            seq = pp.get_sequence()
            chain_sequence = pp.get_sequence()      
            seq_list.append(chain_sequence)
            chain_start_pos = pp[0].get_id()[1] 
            sequences[seq] = chain_start_pos
        for i in seq_list:
            if sequences[i] != 1:
                counter = sequences[i]
            else:
                counter = 0
            for k in i:
                counter +=1
                if k in hydrophobes:
                    for chain in cmd.get_chains(ntpath.basename(self.protein).split('.')[0]):
                        cmd.select('Hydrophobes_' + str(k) + '_'+ str(counter-1) + '_' + str(chain), 'chain ' + str(chain) + ' and  resi  ' + str(counter-1) +  ' and resn ' + str(three_letter[k]) + ' and !name C+N+O')
                        cmd.show('sticks', 'Hydrophobes_' + str(k) + '_'+ str(counter-1) + '_' + str(chain))
                        atom_count = cmd.count_atoms('Hydrophobes_' + str(k) + '_'+ str(counter-1) + '_' + str(chain))
                        if atom_count  > 0 and 'Hydrophobes_' + str(k) + '_'+ str(counter-1) + '_' + str(chain) not in hydrophobe_resi :
                            hydrophobe_resi.append('Hydrophobes_' + str(k) + '_'+ str(counter-1) + '_' + str(chain))
        for resi in hydrophobe_resi:
            print(resi)
            resi_pos = hydrophobe_resi.index(resi)
            for resi2 in hydrophobe_resi:
                dist = cmd.distance('Hydrophobe_contacts_' + str(resi) + '_' + str(resi2), resi, resi2, '4', mode=0)
                if dist >0:
                    D = get_raw_distances('Hydrophobe_contacts_' + str(resi) + '_' + str(resi2))
                    for i in D:
                        hydrophobe_dict.setdefault('Hydrophobe_contacts_' + str(resi) + '_' + str(resi2), []).append(i)
                cmd.delete('Hydrophobe_contacts_' + str(resi) + '_' + str(resi2))
        shortest_bond = {}
        hydro_shrt_bond = {}
        for key in hydrophobe_dict:
            for bond_data in hydrophobe_dict[key]:
                bond_info = list(bond_data)
                shortest_bond.setdefault(key, []).append(bond_info[2])
        for bonds in shortest_bond:
            shortest_bond[bonds].sort()
        for key in shortest_bond:
            shortest_bond[key] = shortest_bond[key][0]
        for key2 in shortest_bond:
            for key3 in hydrophobe_dict:
                for record in hydrophobe_dict[key3]:
                    if shortest_bond[key2] in record:
                        for data in record:
                            if ntpath.basename(self.protein).split('.')[0] in str(data):
                                    hydro_shrt_bond.setdefault(key2, []).append(list(str(data[1])))
        for atom_data in hydro_shrt_bond:
            atom_bond = []
            for data in hydro_shrt_bond[atom_data]:
                atom_pos = ''
                for num in data:
                    atom_pos+= num
                atom_bond.append(atom_pos)
            hydrophobe_dict[atom_data] = atom_bond
        for bond in hydrophobe_dict:
            cmd.select(bond + '_atom1', 'id ' + str(hydrophobe_dict[bond][0]))
            cmd.select(bond + '_atom2', 'id ' + str(hydrophobe_dict[bond][1]))
            cmd.distance('Hydrophobic_bond_' + bond + '_atom1_' +  bond + '_atom2', bond + '_atom1', bond + '_atom2',  '4', mode = 0)
            cmd.color('orange', 'Hydrophobic_bond_' + bond + '_atom1_' +  bond + '_atom2')
        
    def h_bonds(self, model):
        three_letter ={'V': 'VAL', 'I': 'ILE', 'L': 'LEU', 'E': 'GLU', 'Q': 'GLN', 'D': 'ASP', 'N': 'ASN', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE', 'Y': 'TYR', 'R': 'ARG', 'K': 'LYS', 'S': 'SER', 'T': 'THR', 'M': 'MET', 'A': 'ALA', 'G': 'GLY', 'P': 'PRO', 'C': 'CYS'}
        self.protein = model
        prot_path = os.path.dirname(os.path.abspath(self.protein))
        os.chdir(prot_path)
        parser = PDBParser(PERMISSIVE=1)
        structure_id = ntpath.basename(self.protein).split('.')[0]
        filename = ntpath.basename(self.protein)
        structure = parser.get_structure(structure_id, filename)
        ppb=PPBuilder()
        ppb.build_peptides(structure)
        cmd.load(self.protein)
        cmd.show('sticks')
        cmd.hide('cartoon')
        sequences = {}
        for chain in cmd.get_chains(ntpath.basename(self.protein).split('.')[0]):
            cmd.select('Chain_' + str(chain), 'chain ' + str(chain))
            for chain2 in cmd.get_chains(ntpath.basename(self.protein).split('.')[0]):
                cmd.select('Chain2_' + str(chain2), 'chain ' + str(chain2))
                cmd.distance("H-bonds_" + str(chain) + str(chain2), 'Chain_' + str(chain), 'Chain2_' + str(chain2), '4', mode=2)
                cmd.color('cyan', "H-bonds_" + str(chain) + str(chain2))
                
                
            
p = BOND_ID()
p.ionic_bonds('/home/rotan/Desktop/7ntk_D.pdb')
p.hydrophobic_bonds('/home/rotan/Desktop/7ntk_D.pdb')
p.h_bonds('/home/rotan/Desktop/7ntk_D.pdb')