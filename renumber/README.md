```renumber_re.py``` will renumber the index of amino acid in ```.pdb``` to the same index in ```.fasta```(both from RCSB).  
You need to provide the addresses of fasta file and pdb file,name of the chain that need to be renumbered, and the saving address. Read fasta in ```str``` and pdb using ```PDBParser()``` from ```BioPython```.

Example Code Snippet  
````
fasta = next(
        SeqIO.parse('/projectnb2/docking/imhaoyu/my_protein_toolbox/tests/renumber_res/case_2/7MDE.fa', 'fasta'))
fasta_seq = fasta.seq  
parser = PDBParser()
structure = parser.get_structure('7MDE',
                                     '/projectnb2/docking/imhaoyu/my_protein_toolbox/tests/renumber_res/case_2/7MDE.pdb')    
renumber(fasta_seq, structure, chain_name='A', new_PDB_name='7MDE_AB')    
````
We renumber the index of chain A in 7MDE.pdb and save it in 7MDE_AB.pdb .

There are cases that the amino-acid are non-standard (for example, selenomethionines(MSE)). 
