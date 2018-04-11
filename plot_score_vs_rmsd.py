#! /usr/bin/env python2.7
from optparse import OptionParser
import matplotlib
import matplotlib.pyplot as plt
plt.style.use('ggplot')
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import sys

__author__ = 'Alex Sevy'

def normalize_data(x1,y1):
	x,y = zip(*sorted(zip(x1,y1), key=lambda a: a[0]))
	p_95 = np.percentile( y, 95)
	p_5 = np.percentile( y, 5)
	min_value = min(y)
	max_value = max(y)
	y=[float(z - min_value) / (p_95 - p_5) for z in y]
	return x,y

# takes normalized and sorted rms/score lists
def disc( rms1, score1 ):
	rms, score = normalize_data( rms1, score1 )
	rms_bins = [0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.5,4.0,5.0,6.0]
	minimum_bin_population = 20

	rmses = np.array( rms )
	normalized_scores = np.array( score )

	bin_scores = []
	for bin in rms_bins:

		below_bin = [sc for rm,sc in zip(rms, score) if rm < bin]
		above_bin = [sc for rm,sc in zip(rms, score) if rm > bin]

		if len(below_bin) < minimum_bin_population or \
		   len(above_bin) < minimum_bin_population:
			bin_scores.append( 0.0 )
			continue

		min_below_bin = min( below_bin )
		min_above_bin = min( above_bin )

		gap = min_below_bin - min_above_bin
		bin_scores.append( gap )

	return bin_scores

if __name__ == '__main__':
	usage = "%prog [options] <score_vs_rmsd_table>"
	parser=OptionParser(usage)
	parser.add_option("--title",dest="title", help="title of the plot", default="")
	parser.add_option('--modelkey',dest='modelkey', help='only plot models with this key in their name',default='')
	parser.add_option('--color',dest='color', help='marker color for the plot',default='')
	parser.add_option('--no-funnel-score',dest='no_funnel_score', help='should I remove the legend with the funnel score', action='store_true', default=False)
	parser.add_option('--font-size',dest='font_size', help='font size', default=12 )
	parser.add_option('--save',dest='save', help='Should I save this figure?', action='store_true', default=False)
	(options,args)= parser.parse_args()

	rmsd_min = 0
	rmsd_max = 0
	score_min = 0
	score_max = -1000000

	for index, arg in enumerate(args):
		lines = open(arg).readlines()
		
		header = lines[0].split()
		score_index = int(header.index('score'))
		rmsd_index = int(header.index('RMSD'))
		rmsd,score = zip(*[(float(line.split()[rmsd_index]), float(line.split()[score_index])) for line in lines[1:] if options.modelkey in line.split()[-1]])
		rmsd_min = min(rmsd_min,min(rmsd))
		rmsd_max = max(rmsd_max,max(rmsd))
		score_min = min(score_min,min(score))
		score_max = max(score_max,max(score))

	for index, arg in enumerate(args):
		lines = open(arg).readlines()
		
		header = lines[0].split()
		core_index = int(header.index('score'))
		rmsd_index = int(header.index('RMSD'))
		rmsd,score = zip(*[(float(line.split()[rmsd_index]), float(line.split()[score_index])) for line in lines[1:] if options.modelkey in line.split()[-1]])

		plt.rc('font',size=int(options.font_size))
		plt.rc('axes', titlesize=int(options.font_size))
		plt.rc('axes', labelsize=int(options.font_size))
		plt.figure(index)
		
		color = options.color
		if color == '':
			color = 'b' if index == 0 else 'r'
		plt.plot( rmsd, score, 'o', color=color, alpha=0.75 )
		plt.axis((-0.1,rmsd_max+2,score_min-2,score_max+2))
		plt.xlabel( 'CA RMSD' )
		plt.ylabel( 'Score' )
		if options.title:
			title = options.title
		else:
			title = arg.split('/')[-1].split('.')[0]
		plt.title( title )

		if not options.no_funnel_score:
			funnel_discrimination = disc( rmsd, score )

			plt.legend(['Funnel score: %0.2f'%np.mean(funnel_discrimination)], loc='upper right', numpoints=None, fontsize='medium')
	
		## if this is the last plot wait before exiting
		if options.save:
			with PdfPages(arg.split('.')[0]+'.pdf') as pp:
				pp.savefig()
				plt.close()
		else:	
			if index == len(args)-1:
				plt.show()

			## otherwise keep moving
			else:
				plt.show(block=False)

