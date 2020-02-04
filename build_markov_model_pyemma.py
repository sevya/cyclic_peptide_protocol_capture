#! /usr/bin/env python

import os
import numpy as np
import pyemma
from pyemma.util.contexts import settings
import mdtraj as md
from argparse import ArgumentParser
import mdtraj

# Main - used to load data and execute functions
def main():
	# Get commandline parser
	args = cmdlineparse()

	# Run Featurize
	features = Featurize(args)

	# Store features to reader
	reader = pyemma.coordinates.source(args.trajfiles, features=features)

	# Convert trajectory to time-lagged independent component analysis (TICA) feature trajectory
	tica = pyemma.coordinates.tica(reader, lag=int(args.lag), var_cutoff=float(args.tica))
	feat_traj = tica.get_output(stride=int(args.stride))

	# Cluster and estimate maximum likelihood markov model
	cluster = pyemma.coordinates.cluster_kmeans(tica, k=int(args.k_clusters), max_iter=int(args.k_iter))
	msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=int(args.lag), dt_traj=str(args.dt))

	# Check reversible connectivity
	print('fraction of states used = {:f}'.format(msm.active_state_fraction))
	print('fraction of counts used = {:f}'.format(msm.active_count_fraction))

	# Run Plot
	SampleStates(msm, args)

# Generate feature trajectories
def Featurize(ARGS):

	# Load topology file and generate feature object
	features = pyemma.coordinates.featurizer(ARGS.topfile)

	# Add dihedral features
	if ARGS.backbone_torsions:
		features.add_backbone_torsions(deg=True, cossin=True, periodic=ARGS.periodic)
	# Add sidechain features
	if ARGS.sidechain_torsions:
		features.add_sidechain_torsions(deg=True, cossin=True, which='all', periodic=ARGS.periodic)

	# Print feature labels and return
	print(features.describe())
	return features

def SampleStates(MSM,ARGS):
	"Sample from metastable states"
	nstates=int(ARGS.n_metastable_states)
	MSM.pcca(nstates)

	# Print stationary distributions
	ofile = open(str(ARGS.output_prefix)+".stationary_distributions.dat", "w+")
	for i, s in enumerate(MSM.metastable_sets):
		ofile.write('π_{} = {:f}'.format(i, MSM.pi[s].sum()))
		ofile.write('\n')
	ofile.close()

	index = int(0)
	for idist in MSM.sample_by_distributions(MSM.metastable_distributions, int(ARGS.n_samples_each_state)):
		my_samples = pyemma.coordinates.save_traj(ARGS.trajfiles,idist,str(ARGS.output_prefix)+"."+str(index)+".metastable.nc",top=str(ARGS.topfile),image_molecules=bool(ARGS.periodic))
		index = index + 1

# Argument parser
def cmdlineparse():
	# Construct argument parser
	parser = ArgumentParser(description="Command line arguments",add_help=True)

	# Make multiple argument groups for organization purposes
	global_group = parser.add_argument_group("Global","General arguments")
	feature_group = parser.add_argument_group("Features","Indicate what types of features to generate for the trajectory")
	dimensions_group = parser.add_argument_group("Dimensionality Reduction and Clustering", "Arguments related to dimensionality reduction and clustering methods")

	# Global group
	global_group.add_argument("-topfile", dest="topfile", required=True, help="Topology file", metavar="<topology file>")
	global_group.add_argument("-trajfiles", dest="trajfiles", nargs="+", required=False, help="Trajectory files",metavar="<trajectory files>")
	global_group.add_argument("-output_prefix", dest="output_prefix", required=False,default="PyemmaOutput", help="Output prefix for output files",metavar="<output files>")
	global_group.add_argument("-stride", dest="stride", required=False, default=1, help="Load every <stride>th trajectory frame", metavar="<stride>")
	global_group.add_argument("-dt", dest="dt", required=False, default='1 step', help="Frame step quantity and unit", metavar="<dt>")
	global_group.add_argument("-lag", dest="lag", required=False, default='1', help="Number of lag steps", metavar="<lag>")
	global_group.add_argument("-periodic",dest="periodic",required=False, default=False, action="store_true", help="Indicate the presence of periodic boundary conditions")

	# Feature group
	feature_group.add_argument("-backbone_torsions",dest="backbone_torsions",required=False, action="store_true", default=False,
						help="Add all backbone torsions to the feature list")
	feature_group.add_argument("-sidechain_torsions",dest="sidechain_torsions",required=False, action="store_true", default=False,
						help="Add all sidechain torsions to the feature list")

	# Dimensions group
	dimensions_group.add_argument("-tica",dest="tica",required=False,
								  help="Perform TICA with the requested kinetic variance cutoff",
								  metavar="<tica_cutoff>")
	dimensions_group.add_argument("-k_clusters",dest="k_clusters",required=False,default=100,
								  help="Number of clusters for k-means clustering",
								  metavar="<k_klusters>")
	dimensions_group.add_argument("-k_iter",dest="k_iter",required=False,default=10,
								  help="Number of iterations for k-means clustering",
								  metavar="<k_iter>")
	dimensions_group.add_argument("-n_metastable_states", dest="n_metastable_states", required=False, help="Number of metastable states for network analysis", metavar="<n_metastable_states>")
	dimensions_group.add_argument("-n_samples_each_state", dest="n_samples_each_state", required=False, help="Number of structures to sample from each metastable state", metavar="<n_samples_each_state>")

	# Done
	args=parser.parse_args()
	return args

# Run from commandline
if __name__ == '__main__':
	main()

