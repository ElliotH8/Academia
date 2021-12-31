def mass(file1):
    file = open(file1, 'r')
    lines = file.readlines()
    file.close()
    
    mass_dict = {'ALA': 59, 'ARG': 89, 'ASN': 115, 'ASP': 122, 'CYS': 155, 'GLU': 117, 'GLN': 117, 'GLY': 105, 'HIS': 109, 'ILE': 130, 'LEU': 50, 'LYS': 45, 'MET': 43, 'PHE': 42, 'PRO': 51, 'SER': 27, 'THR': 32, 'TRP': 105, 'TYR': 99, 'VAL': 92}
    
    aatype = []
    #type of atom
    amino_acid = []
    #amino acid atom belongs to
    total_aminoacids = []
    #list of all amino acids
    
    for element in lines:
        aatype.append(element.split()[2])
        amino_acid.append(element.split()[3])
    
    for iteration, element in enumerate(aatype):
        if element == 'CA':
            m = amino_acid[iteration]
            total_aminoacids.append(m)
        else:
            continue
        
    protein_mass = []
    #total mass of protein
    for element in total_aminoacids:
        protein_mass.append(mass_dict[element]) 

    return print(sum(protein_mass))
#protein string function

def proteinstring(file1):
    file = open(file1, 'r')
    lines = file.readlines()
    file.close()
    
    short_dict = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
    aatype = []
    #type of atom
    amino_acid = []
    #amino acid atom belongs to
    total_aminoacids = []
    #list of all amino acids
    
    for element in lines:
        aatype.append(element.split()[2])
        amino_acid.append(element.split()[3])
    
    for iteration, element in enumerate(aatype):
        if element == 'CA':
            m = amino_acid[iteration]
            total_aminoacids.append(m)
        else:
            continue
        
    protein_name = []
    #list of amino acids converted to one letter
    for element in total_aminoacids:
        protein_name.append(short_dict[element])
    protein_str = ''
    for element in protein_name:
        protein_str += str(element)
    return print(protein_str)

mass('1ubq.pdb')
proteinstring('1ubq.pdb')
