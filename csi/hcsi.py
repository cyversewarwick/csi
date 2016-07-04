import numpy as np
import scipy as sp
import multiprocessing as mp
import pandas as pd
import argparse
import sys
import warnings
import pickle

import csi
import AbstractGibbs_hCSI as ag
from main import parse_gp_hyperparam_priors

def loadData(fd):
	'''
	Import a CSI-formatted CSV file. Based on Sam's loadData() from his CSI implementation.
	
	Input:
		* fd - argparse-generated file handle for the CSV file
	'''
	# read in data and make sure headers are correct
	inp = pd.read_csv(fd,dtype=str,index_col=0,header=[0,1])
	inp.columns = pd.MultiIndex.from_tuples([(a,float(b)) for a,b in inp.columns],names=inp.columns.names)
	# convert to floating point values
	return inp.astype(float)

def parse_args():
	'''
	Parse the command line arguments provided and return an args structure with them inside.
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument('--Input', dest='input', type=argparse.FileType('r'), required=True, help='The CSV file featuring your expression data. First column for gene names, first row for condition names repeated for each condition time point, second row for condition time point information.')
	parser.add_argument('--Depth', dest='depth', type=int, default=2, help='CSI parental set depth truncation to maintain computational tractability. Default: 2')
	parser.add_argument('--Prior', dest='gpprior', type=str, default='10;0.1', help='CSI Gaussian Process prior, provided as \'shape;scale\' for a gamma distribution or \'uniform\' for a uniform distribution. Default: \'10;0.1\'')
	parser.add_argument('--BetaPrior', dest='betaprior', type=str, default='1;1', help='hCSI temperature prior, provided as \'shape;scale\' for a gamma distribution. Default: \'1;1\'')
	parser.add_argument('--SeedOffset', dest='offset', type=int, default=0, help='Shift the seeding (based on child gene row number in the CSV) by a fixed amount to get (slightly) different run results. Default: 0 (no offset)')
	parser.add_argument('--Normalise', dest='normalise', choices=['none','center','standardise'], default='standardise', help='Choice of normalisation - \'none\' to leave data as is, \'center\' to center, \'standardise\' to standardise. Default: \'standardise\'')
	parser.add_argument('--Genes', dest='genes', default=None, nargs='+', help='Child gene set to evaluate, if you wish to only run hCSI on a subset of the available gene space. Provide as space delimited names matching the CSV file. Default: None (analyse the whole dataset)')
	parser.add_argument('--Pool', dest='pool', type=int, default=0, help='Number of threads to open up for parallelising hCSI on a per-gene basis. Default: 0 (automated parallelising)')
	parser.add_argument('--LikelihoodPool', dest='likpool', default=None, type=int, help='Likelihood computation is the most resource-intensive part of hCSI. This option allows an additional level of parallelisation in likelihood computation, to speed up run times on resource-heavy local computational platforms. Default: None (no additional parallelisation level)')
	parser.add_argument('--Samples', dest='samples', type=int, default=25000, help='Number of Gibbs updates to perform within hCSI. Default: 25,000')
	parser.add_argument('--BurnIn', dest='burnin', type=int, default=2500, help='Number of initial Gibbs updates to discard as burn-in. Default: 2,500')
	parser.add_argument('--Pickle', dest='pickle', action='store_true', help='Flag. If provided, the obtained Gibbs value chains for individual models and the hypernetwork are stored as a Python Pickle. Refer to readme for more in depth formatting information.')
	parser.add_argument('--GibbsVerbose', dest='verbose', action='store_true', help='Flag. If provided, the script will print status updates once every 5,000 Gibbs updates.')
	args = parser.parse_args()
	return args

def hamming(p1, p2):
	'''
	Compute the Hamming distance between two parental sets.
	
	Input:
		* p1, p2 - the two parental sets, standard tuple formatting, i.e. ([parents], child)
	'''
	#filter out to what we actually want from the parents
	p1_set = set(p1[0])
	p2_set = set(p2[0])
	#the hamming distance is the total of elements that differ between the parent combos
	ulen = len(p1_set.union(p2_set))
	ilen = len(p1_set.intersection(p2_set))
	#so we return the total number that only show up in one of the two
	return ulen-ilen

def ismember(a, b):
	'''
	Functional clone of Matlab's ismember(), returning a Boolean mask indicating which elements of a are in b.
	
	Input:
		* a, b - the two lists to be compared. Need to be lists, even if a or b are single elements.
	'''
	bind = {}
	for i, elt in enumerate(b):
		if elt not in bind:
			bind[elt] = True
	return np.asarray([bind.get(itm, False) for itm in a])

def roulette(dist):
	'''
	Roulette wheel sampling algorithm from a distribution. Returns the index of the sampled position.
	
	Input:
		* dist - the pdf of the distribution to be sampled from, ensure it adds up to 1
	'''
	dist2 = np.cumsum(dist)
	return np.argmax(dist2 >= sp.rand(1))

class GibbsHCSI(ag.AbstractGibbs):
	def __init__(self, rvObjList, rvHyperNetwork):
	#add a separate hypernetwork rv
		self.rvl = rvObjList
		self.indexList = range(len(self.rvl))
		self.sampledValues = []
		for i in self.indexList:
			self.sampledValues.append([])
		self.rvHyperNetwork = rvHyperNetwork
		self.sampledValuesHyperNetwork = []

	def gibbsUpdate(self,iteration_number):
	#expand the thing with hypernetwork sampling and parameter MCMC
		#run through all the variables
		valList = []
		betaList = []
		pset = self.rvHyperNetwork.getValRange()
		for i in self.indexList:
			#for each variable, build a distribution accross all possible values for that variable
			#based on current values of other variable 
			distribution = self.rvl[i].getConditionalDistribution(self.rvHyperNetwork.getCurrentValue())
			#sample a random value based on the distribution
			valInd = roulette(distribution)
			self.rvl[i].setCurrentValue(pset[valInd])
			valList.append(pset[valInd])
			betaList.append(self.rvl[i].beta)
		#sample the hypernetwork
		distribution = self.rvHyperNetwork.getConditionalDistribution(valList,betaList)
		valInd = roulette(distribution)
		self.rvHyperNetwork.setCurrentValue(pset[valInd])
		#hyperparameter and temperature sampling
		for i in self.indexList:
			#sample the hyperparameters
			old_hypers = self.rvl[i].em.hypers
			new_hypers = sp.exp(sp.log(old_hypers)+self.rvl[i].thetaconst*sp.randn(3))
			#set up the thetajumps (nohypers to not shift RNG chain) and compute logliks
			self.rvl[i].thetajump.setup([valList[i]])
			self.rvl[i].thetajump.hypers = old_hypers
			old_ll = self.rvl[i].thetajump.logliks()[0]
			self.rvl[i].thetajump.hypers = new_hypers
			new_ll = self.rvl[i].thetajump.logliks()[0]
			#numeric thing before exping, scale to largest likelihood
			(new_ll,old_ll) = (new_ll,old_ll)-np.max((new_ll,old_ll))
			#compute P(hypers)
			p_hypers_old = np.prod(sp.stats.gamma.pdf(old_hypers,a=self.rvl[i].em._prior_shape,scale=self.rvl[i].em._prior_scale))
			p_hypers_new = np.prod(sp.stats.gamma.pdf(new_hypers,a=self.rvl[i].em._prior_shape,scale=self.rvl[i].em._prior_scale))
			#transition probability
			p_accept_theta = (sp.exp(new_ll) * p_hypers_new) / (sp.exp(old_ll) * p_hypers_old)
			if np.random.uniform() <= p_accept_theta:
				#we accept the new hypers
				self.rvl[i].em.hypers = new_hypers
				self.rvl[i].thetacount += 1
			'''VISIBLE BREAK'''
			#sample the beta (apparently Chris doesn't use the log space shift thing here)
			old_beta = self.rvl[i].beta
			new_beta = old_beta+self.rvl[i].betaconst*sp.randn()
			#prepare the scaling constants
			Z = np.zeros(len(self.rvl[i].em.pset))
			for j in np.arange(len(self.rvl[i].em.pset)):
				Z[j] = hamming(self.rvl[i].em.pset[j],self.rvHyperNetwork.currentValue)
			const_old = np.sum(np.exp((-1)*old_beta*Z))
			const_new = np.sum(np.exp((-1)*new_beta*Z))
			#compute the probabilities
			p_beta_old = sp.exp((-1)*old_beta*hamming(valList[i],self.rvHyperNetwork.currentValue))/const_old * sp.stats.gamma.pdf(old_beta,a=self.rvl[i].betaprior[0],scale=self.rvl[i].betaprior[1])
			p_beta_new = sp.exp((-1)*new_beta*hamming(valList[i],self.rvHyperNetwork.currentValue))/const_new * sp.stats.gamma.pdf(new_beta,a=self.rvl[i].betaprior[0],scale=self.rvl[i].betaprior[1])
			p_accept_beta = p_beta_new/p_beta_old
			if np.random.uniform() <= p_accept_beta:
				#we accept the new beta
				self.rvl[i].beta = new_beta
				self.rvl[i].betacount += 1
			'''VISIBLE BREAK'''
			#re-evaluate constants on 100 Gibbs goes
			if iteration_number % 100 == 0:
				#aim for 25 jumps on each parameter because stats and Chris said so
				#too many jumps, too local, push it out of its comfort zone
				if self.rvl[i].thetacount > 35:
					self.rvl[i].thetaconst *= 1.1
				#too few jumps, hopping around the place too much, localise it a bit
				elif self.rvl[i].thetacount < 15:
					self.rvl[i].thetaconst *= 0.9
				if self.rvl[i].betacount > 35:
					self.rvl[i].betaconst *= 1.1
				elif self.rvl[i].betacount < 15:
					self.rvl[i].betaconst *= 0.9
				self.rvl[i].betacount = 0
				self.rvl[i].thetacount = 0

	def sample(self, repeats=1, verbose=False):
	#add hypernetwork sample storage
		for i in range(repeats):
			self.gibbsUpdate(i+1)
			for j in self.indexList:
				self.sampledValues[j].append(self.rvl[j].getCurrentValue())
			self.sampledValuesHyperNetwork.append(self.rvHyperNetwork.getCurrentValue())
			if (i+1) % 5000 == 0:
				if verbose:
					print(str(i+1)+' iterations done.')

class RandomVariableCondition(ag.RandomVariable):
	def __init__(self, csidata, cond, gene, gpprior, betaprior, depth, standardise):
	#override init to include all sorts of CSI stuffs
		#extract the data frame columns that actually have our condition's data
		#(without killing the fine balance of the tuple indexing)
		colNames = np.asarray([x[0] for x in csidata.columns.values])
		numtemp = np.arange(len(colNames))
		inp = csidata.iloc[:,numtemp[colNames==cond]]
		#standardising
		if standardise=='standardise':
			inp[:][:] = sp.stats.mstats.zscore(inp,axis=1,ddof=1)
		elif standardise=='center':
			inp[:][:] = inp[:][:] - np.mean(inp[:][:],axis=1)[:,None]
		#some processing stuff borrowed from the CSI code
		assert (inp.columns.is_monotonic_increasing)
		self.cc = csi.Csi(inp)
		self.em = self.cc.getEm()
		#override exp(N(0,1)) hypers with U(0,1) hypers
		self.em.hypers = sp.rand(3)
		if gpprior:
			self.em.set_priors(gpprior[0], gpprior[1])
		#prepare the EM object
		self.em.sampleinitweights = False
		self.em.setup(self.cc.allParents(gene,depth),args.likpool)
		#beta initialised at 0.1 as per Chris code (1 was his recommendation)
		self.beta = 0.1
		self.betaprior = betaprior
		#MCMC constants (0.1 in code, 1 recommendation)
		self.thetaconst = 0.1
		self.betaconst = 0.1
		self.thetacount = 0
		self.betacount = 0
		#random initialisation
		ind = roulette(np.ones(len(self.em.pset))/len(self.em.pset))
		self.currentValue = self.em.pset[ind]
		self.valRange = self.em.pset
		self.distribution = list(range(len(self.valRange)))
		#thetajump setup, everything will be dealt with later
		self.thetajump = self.cc.getEm()
		self.thetajump.sampleinitweights = False

	def getConditionalDistribution(self, hyperparent):
	#override with actual computation
		i = 0
		logliks = self.em.logliks()
		#turn logliks into actual numbers, scaled to the largest one
		logliks = np.exp(logliks-np.max(logliks))
		#for each possible value of the RV
		for v in self.getValRange():
			#skipping P(theta) and all the Z's as they're the same for each parent combo
			p =  logliks[i] * np.exp((-1)*self.beta*hamming(v,hyperparent))
			#set the distribution
			self.distribution[i] = p
			i = i + 1
		self.normalizeDistribution()
		return self.distribution

	def normalizeDistribution(self):
	#we actually need the normalisation
		self.distribution = self.distribution/np.sum(self.distribution)

class RandomVariableHyperNetwork(ag.RandomVariable):
	def __init__(self, csidata, gene, depth):
	#override init to include all sorts of CSI stuffs
		cc = csi.Csi(csidata)
		self.valRange = cc.allParents(gene,depth)
		ind = roulette(np.ones(len(self.valRange))/len(self.valRange))
		self.currentValue = self.valRange[ind]
		self.distribution = list(range(len(self.valRange)))

	def getConditionalDistribution(self, condParents, betas):
	#override with actual computation
		i = 0
		for v in self.getValRange():
			p = 1
			for (j, cond) in enumerate(condParents):
				p = p * np.exp((-1)*betas[j]*hamming(v,cond))
			self.distribution[i] = p
			i = i + 1
		self.normalizeDistribution()
		return self.distribution

	def normalizeDistribution(self):
	#we actually need the normalisation
		self.distribution = self.distribution/np.sum(self.distribution)

def runGibbs(gene_id, inp, gpprior, betaprior, args):
	'''
	Perform a complete hCSI run for a given child gene, returning the relevant rows of the marginal matrix.
	
	Input:
		* gene_id - which row of the input CSV is the child to be analysed, used for seed setting and data access
		* inp - the parsed CSV expression file as a Pandas data frame
		* gpprior - [shape, scale] information on the gamma prior for the CSI Gaussian process hyperparameters
		* betaprior - [shape, scale] information on the gamma prior for the temperature parameters
		* args - the parsed command line arguments
	'''
	#fish out the single gene via gene_id, also set the seeds via that
	gene = args.genes[gene_id]
	np.random.seed(gene_id+args.offset)
	#commence proper part of thing
	print('Processing '+gene+'...')
	hnrv = RandomVariableHyperNetwork(inp,gene,args.depth)
	crv = []
	conditions = np.unique([x[0] for x in inp.columns.values])
	for cond in conditions:
		crv.append(RandomVariableCondition(inp,cond,gene,gpprior,betaprior,args.depth,args.normalise))
	gibbs = GibbsHCSI(crv,hnrv)
	gibbs.sample(args.samples, args.verbose)
	#ditch the burn-in
	for i in np.arange(len(gibbs.sampledValues)):
		gibbs.sampledValues[i] = gibbs.sampledValues[i][args.burnin:]
	gibbs.sampledValuesHyperNetwork = gibbs.sampledValuesHyperNetwork[args.burnin:]
	#parse rows of the magical marginal matrices
	out = np.zeros((len(conditions)+1,len(args.genes)))
	numtemp = np.arange(len(args.genes))
	for i in np.arange(len(conditions)):
		for j in np.arange(len(gibbs.sampledValues[i])):
			mask = ismember(args.genes, list(gibbs.sampledValues[i][j][0]))
			out[i,numtemp[mask]] += 1
		out[i,:] /= len(gibbs.sampledValues[i])
	#repeat for the hypernetwork
	for i in np.arange(len(gibbs.sampledValuesHyperNetwork)):
		mask = ismember(args.genes, list(gibbs.sampledValuesHyperNetwork[i][0]))
		out[-1,numtemp[mask]] += 1
	out[-1,:] /= len(gibbs.sampledValuesHyperNetwork)
	#potential chain storing
	if args.pickle:
		chain = gibbs.sampledValues
		chain.append(gibbs.sampledValuesHyperNetwork)
	else:
		chain = None
	return (out, chain, gene)

def _pool_init(inp_, gpprior_, betaprior_, args_):
	'''
	Helper function to set up the "global" variables for the per-gene parallelised pool
	'''
	global inp, gpprior, betaprior, args
	inp = inp_
	gpprior = gpprior_
	betaprior = betaprior_
	args = args_

def pool_runGibbs(gene_id):
	'''
	Helper function calling runGibbs() using the "global" variables made by _pool_init
	
	Input:
		* gene_id - which row of the input CSV is the current child, the only variable part of the input
	'''
	return runGibbs(gene_id, inp, gpprior, betaprior, args)

def main():
	#read arguments
	args = parse_args()
	#parse data
	inp = loadData(args.input)
	conditions = np.unique([x[0] for x in inp.columns.values])
	#assess indegree sanity
	if args.depth < 1:
		sys.stderr.write("Error: truncation depth must be greater than or equal to one")
		sys.exit(1)
	#assess parallelisation
	if args.pool == 0:
		args.pool = mp.cpu_count()
	#parse priors
	if args.gpprior is None or args.gpprior == 'uniform':
		gpprior = None
	else:
		try:
			gpprior = parse_gp_hyperparam_priors(args.gpprior)
		except ValueError(s):
			sys.stderr.write("Error: "+s)
			sys.exit(1)
	try:
		betaprior = parse_gp_hyperparam_priors(args.betaprior)
	except ValueError(s):
		sys.stderr.write("Error: "+s)
		sys.exit(1)
	#parse gene list
	genes = args.genes
	if genes is None:
		genes = list(inp.index)
		#not pretty, but passes it on to the parallelised workers with ease
		args.genes = genes
	else:
		missing = np.setdiff1d(genes, inp.index)
		if len(missing) > 0:
			sys.stderr.write("Error: The following genes were not found: {missing}\n".format(missing=', '.join(missing)))
			sys.exit(1)
	#prepare output things
	output = []
	for i in range(len(conditions)+1):
		template = np.ndarray((len(genes)+1,len(genes)+1),dtype=object)
		template[0,0] = ''
		template[1:,0] = np.asarray(genes)
		template[0,1:] = np.asarray(genes)
		output.append(template)
	numtemp = 1+np.arange(len(genes))
	chains = []
	#"Parallelised or not, here I come!" - Gibbs sampler, 2016
	if args.pool > 1:
		p = mp.Pool(args.pool, _pool_init, (inp, gpprior, betaprior, args))
		for (out, chain, gene) in p.imap_unordered(pool_runGibbs, np.arange(len(genes))):
			#wrap it into a one-element list so that ismember sees it whole
			mask = ismember(genes,[gene])
			ind = numtemp[mask][0]
			#the code will spit out Nones if not set to pickle. append with peace
			chains.append(chain)
			#sneaking in the hypernetwork as the last "condition"
			for i in np.arange(len(conditions)+1):
				output[i][ind,1:] = np.asarray([str(x) for x in out[i,:]])
	else:
		for gene_id in genes:
			(out, chain, gene) = runGibbs(gene_id, inp, gpprior, betaprior, args)
			#wrap it into a one-element list so that ismember sees it whole
			mask = ismember(genes,[gene])
			ind = numtemp[mask][0]
			#the code will spit out Nones if not set to pickle. append with peace
			chains.append(chain)
			#sneaking in the hypernetwork as the last "condition"
			for i in np.arange(len(conditions)+1):
				output[i][ind,1:] = np.asarray([str(x) for x in out[i,:]])
	#spit the things out
	if args.pickle:
		dumpfile = open('chains.pickle','wb')
		pickle.dump(chains,dumpfile)
		dumpfile.close()
	conditions = list(conditions)
	conditions.append('hypernetwork')
	for i in np.arange(len(conditions)):
		writer = open('hcsi_'+conditions[i]+'.csv','w')
		for j in np.arange(output[i].shape[0]):
			writer.write(','.join(output[i][j,:])+'\n')
		writer.close()

if __name__ == '__main__':
    main()