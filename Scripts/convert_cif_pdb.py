import os
from Bio.PDB import MMCIFParser, PDBIO

# Set the directory where your .cif files are located
cif_directory = './'  # Update this path to where your files are

# Initialize the parser and the PDBIO writer
parser = MMCIFParser(QUIET=True)
io = PDBIO()

# Loop through all files in the directory
for file_name in os.listdir(cif_directory):
    if file_name.endswith(".cif"):
        cif_file_path = os.path.join(cif_directory, file_name)
        
        # Parse the CIF file
        structure = parser.get_structure(file_name, cif_file_path)
        
        # Set the output file name
        pdb_file_name = file_name.replace(".cif", ".pdb")
        pdb_file_path = os.path.join(cif_directory, pdb_file_name)
        
        # Save the structure as a PDB file
        io.set_structure(structure)
        io.save(pdb_file_path)

        print(f"Converted {file_name} to {pdb_file_name}")
