#!/usr/bin/env python2.7

from Bio import SeqIO
from sevy_utils import Utils
import multiprocessing as mp
from optparse import OptionParser 

__author__ = 'Alex Sevy'

if __name__ == '__main__':
	usage = "%prog [options] <fasta_file>"
	parser=OptionParser(usage)

	parser.add_option('--thread-start',dest='start',help='The PDB res ID where you want to start threading',default=1)
	parser.add_option('--chain',dest='chain', help='Chain that you want to thread over', default='A')
	parser.add_option('--multiproc','-m',dest='multiproc',help='Should I run this job on multiple processors?',default=False, action='store_true')
	
	(options,args)= parser.parse_args()

	resfile_header = 'NATAA\nUSE_INPUT_SC\nstart\n'
	
	def write_fasta( record ):
		with open( str(record.id)+'.resfile', 'w' ) as out:
			out.write('## '+str(record.id)+'\n')
			out.write( resfile_header )
			resno = int(options.start)
			for aa in str(record.seq):
				out.write( ' '.join([str(resno), options.chain, 'PIKAA', aa, '\n']))
				resno += 1


	if options.multiproc:
		## Use all available CPUs
		pool = mp.Pool( processes=mp.cpu_count() )
		records = [record for record in SeqIO.parse(open( args[0] ), 'fasta')]
		pool.map( write_fasta, records )
	else:
		for record in SeqIO.parse(open( args[0] ), 'fasta'):
			write_fasta( record )
	

