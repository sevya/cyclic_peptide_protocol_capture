#! /usr/bin/env python2.7

from Bio.PDB import PDBParser
import gzip, sys

__author__ = 'Alex Sevy'

class Utils:

	amino_acids = {
        'ALA': 'A',
        'CYS': 'C',
        'ASP': 'D',
        'GLU': 'E',
        'GLY': 'G',
        'HIS': 'H',
        'ILE': 'I',
        'LYS': 'K',
        'LEU': 'L',
        'MET': 'M',
        'ASN': 'N',
        'PHE': 'F',
        'PRO': 'P',
        'GLN': 'Q',
        'ARG': 'R',
        'SER': 'S',
        'THR': 'T',
        'VAL': 'V',
        'TRP': 'W',
        'TYR': 'Y',
		'TYS': 'Y',
#		'B53': 'U'
	}

	'''
	Returns the one letter code from three letter code 
	'''
	@staticmethod
	def threetoone( aa ):
		return amino_acids[ aa ]

	'''
	Returns the three letter code from one letter code 
	'''
	@staticmethod
	def onetothree( aa ):
		inverted_dict = {value : key for key, value in Utils.amino_acids.items()}
		return inverted_dict[ aa ]

	'''
	Gets BioPDB Structure object from a filename
	takes care of extension handling
	'''
	@staticmethod
	def get_structure( filename ):
		pdbparser = PDBParser( PERMISSIVE=1 )
		ext = filename.split('.')[-1]
		if ext == 'gz':
			handle = gzip.open(filename)
		else:
			handle = open(filename)
		try:
			structure = pdbparser.get_structure( 'X', handle )
			return structure
		except:
			print sys.exc_info()[1]
			print 'Unable to read structure',filename
			exit()

	'''
	Parses a resfile to return an array of designable residues
	Designable residues are in the format [(1, 'A'), (2, 'A'), ...]
	'''
	@staticmethod
	def parse_resfile( filename ):
		designable_residues = []
		start = False

		for line in open( filename ).readlines():
			line = line.strip()
			
			if len(line) == 0:
				continue
			## if the line is commented out, ignore it
			if line[0] == '#':
				continue

			try:
				## If i've already seen the start line,
				## then try to parse the residues
				if start:
					## if there's a residue range specified,
					## use different logic for parsing
					if '-' in line.split():
						res_start, dash, res_end, chain, action = line.split()[:5]
						if 'ALLAA' in action or action == 'NOTAA':
							for resno in range(int(res_start), int(res_end)+1):
								designable_residues.append((int(resno), chain))
					else:
						resno, chain, action = line.split()[:3]
						if 'ALLAA' in action or action == 'NOTAA':
							designable_residues.append((int(resno), chain))
						## if it says PIKAA but there's more than one to choose from
						elif 'PIKAA' in action and len(line.split()[3]) > 1:
							designable_residues.append((int(resno), chain))
				elif line.split()[0] == 'start':
					start = True

			except IndexError: # if there are no items in this line
				continue

		if not designable_residues:
			print 'Error: resfile did not parse correctly. No residues are set as designable'
			exit()	 

		return designable_residues

	'''
	From a Bio.PDB Structure object, gives you the one letter code of the amino acid at
	position 'resno' on chain 'chain'. This method supports use of wildcards - so, if you 
	pass in resno *, it will return a string of the one letter codes of all amino acids on
	the specified chain. Passing wildcard for both chain and resno will return the entire
	protein sequence
	'''
	@staticmethod
	def residue_from_pdb( structure, chain, resno ):
		model = structure[0]
		if chain == '*':
			all_chains = [c.id for c in structure.get_chains()]
			return ''.join( [Utils.residue_from_pdb( structure, c, resno ) for c in all_chains ] )
		if resno == '*':
			name3s = [res.resname for res in model[ chain ].get_residues()]
			name1s = [ Utils.amino_acids[ n ] for n in name3s ]
			return ''.join( name1s )
		name3 = model[ chain ][ resno ].resname
		
		if name3 in Utils.amino_acids:
			return Utils.amino_acids[ name3 ]
		else:
			return 'U'

	'''
	'''
	@staticmethod
	def resno_from_pdb( structure, chain, resno ):
		model = structure[0]
		if chain == '*':
			all_chains = [c.id for c in structure.get_chains()]
			return [Utils.resno_from_pdb( structure, c, resno ) for c in all_chains] 
		if resno == '*':
			return [(res.id[1], chain) for res in model[ chain ].get_residues()]
		
		return [(model[ chain ][ resno ].id[1], chain)]



	'''
	Takes in a PDB file and return the sequence
	Uses a residue list to decide which residues to include in 
	this sequence. Residue list is an array of tuples holding
	the residue number and the chain ID, i.e. [(1, 'A'), (2, 'A'), ...]
	Returns an array of strings
	'''
	@staticmethod
	def sequence_from_pdb( file, residue_list ):
		structure = Utils.get_structure( file )
		current_seq = ''
		model = structure[0]
		for resno, chain in residue_list:
			current_seq += Utils.residue_from_pdb( structure, chain, resno )
		return current_seq

	'''
	Takes in a PDB file and a list of designable residue with possible wildcards
	and returns a list of designable residues in the format [(1, 'A'), ...]

	Returns an array of tuples
	'''
	@staticmethod
	def parse_designable_residues( file, residue_list ):
		structure = Utils.get_structure( file )
		designable_residues = []
		model = structure[0]

		for resno, chain in residue_list:
			for item in Utils.resno_from_pdb( structure, chain, resno ):
				designable_residues.append(item)
		## Collapse array
		return [item for sublist in designable_residues for item in sublist]
		# return designable_residues

	'''
	Takes in an array of tuples, holding sequences and IDs, 
	in the format [(ID, sequence), (ID2, sequence2)...]
	Returns a string of valid FASTA format that can be output 
	directly to a file
	'''
	@staticmethod 
	def fasta_from_sequences( sequences ):
		return '\n'.join( ['>'+name+'\n'+sequence for name, sequence in sequences] )


	'''
	Returns true if this number is within tol of the nearest integer
	'''
	@staticmethod
	def isint( n, tol=1e-6 ):
		closest_int = round( n, 0 )
		return abs( n - closest_int ) < tol










