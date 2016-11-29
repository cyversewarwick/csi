from __future__ import absolute_import, division, print_function

import sys
import optparse
import csv
import re
import itertools as it

import h5py as h5

import logging

import numpy as np
import scipy as sp
import pandas as pd

import multiprocessing as mp

import csi

logger = logging.getLogger('CSI')

def cmdparser(args):
    # create parser objects
    op = optparse.OptionParser()
    out = optparse.OptionGroup(op, "Output Formats")
    compat = optparse.OptionGroup(op, "Compatibility Options")

    op.set_usage("usage: %prog [options] FILE.csv")
    op.set_defaults(
        verbose=False,
        normalise='standardise',
        depth=2,
        numprocs=0,
        depgenes=[],
        tflist=None,
        gpprior='10;0.1')

    # define general parameters
    op.add_option('-v','--verbose',dest='verbose',action='count',
                  help="Increase verbosity (specify twice for more detail)")
    op.add_option('-d', '--depth',dest='depth',type='int',metavar='D',
                  help="Truncate parental set at depth D")
    op.add_option('--gene',dest='genes',action='append',metavar='GENE',
                  help="Limit analysis to a specific gene (repeat for more than one)")
    op.add_option('--mp',dest='numprocs',type='int',
                  help="Number of CPU cores (worker processes) to use")
    op.add_option('--gpprior',dest='gpprior', metavar='PRIORS',
                  help="Gaussian Process Hyperparameters, 'shape;scale'")
    op.add_option('--normalise',dest='normalise', type='choice',
                  choices=['none','center','standardise'],
                  help="How to normalise the data on input ('none', 'center', 'standardise')")
    op.add_option('--tflist',dest='tflist', metavar='FILE',
                  help="A list of TF IDs, one line per TF, to constrain the parental sets to")

    # define output parameters
    op.add_option_group(out)
    out.add_option('--csv',dest='csvoutput',
                   help="write CSV output to FILE", metavar='FILE')
    out.add_option('--hdf5',dest='hdf5output',
                   help="write HDF5 output to FILE", metavar='FILE')

    # define compatibility options
    op.add_option_group(compat)
    compat.add_option('--weighttrunc',type='float', metavar='V',
                      help="Don't evaluate likelihoods whose weight is less than V")
    compat.add_option('--initweights',type='choice',metavar='TYPE',
                      choices=['uniform','weighted'],
                      help="Initialise weights as either 'uniform' or 'weighted'")

    # parse our command line arguments
    return op.parse_args(args)

def parse_gp_hyperparam_priors(s):
    """parse a string of the form 'a,b' or 'a1 a2 a3,b1 b2 b3' and return as
    None or a tuple of numpy arrays"""
    # define some regular expressions to parse the strings
    f = r'[0-9]*\.?[0-9]*(?:e[+-]?[0-9]+)?'
    r1 = re.compile(r'^ *({f}) *$'.format(f=f))
    r3 = re.compile(r'^ *({f}) +({f}) +({f}) *$'.format(f=f))
    # take a string and return a numpy array or raise an exception
    def parse(a):
        m = r1.match(a)
        if m is not None:
            return np.array(m.groups()[0],dtype=float)
        m = r3.match(a)
        if m is not None:
            return np.array(m.groups(),dtype=float)
        raise ValueError("Unable to parse {0} as hyperparameter prior".format(repr(a)))
    # split the string into two
    arr = s.split(';')
    if len(arr) != 2:
        raise ValueError('Expecting only one semicolon in hyperparameters')
    # parse the two sides
    return (parse(arr[0]),
            parse(arr[1]))

def main(args=None):
    if args is None:
        args = sys.argv[1:]

    # parse command line arguments
    (op,fname) = cmdparser(args)

    # extract out the logging level early
    log_level = logging.WARNING
    if   op.verbose == 1: log_level = logging.INFO
    elif op.verbose >= 2: log_level = logging.DEBUG

    # configure logger
    logging.basicConfig(level=log_level)
    logging.getLogger('GP').setLevel(logging.WARNING)
    logging.getLogger('parameters changed meta').setLevel(logging.WARNING)

    # make sure we're only processing a single file
    if len(fname) != 1:
        if len(fname) == 0:
            sys.stderr.write("Error: Please specify the filename to process, or run with '-h' for more options\n")
        else:
            sys.stderr.write("Error: Only one input filename currently supported\n")
        sys.exit(1)

    # pull out the parental set trunction depth and validate
    depth = op.depth
    if depth < 1:
        sys.stderr.write("Error: truncation depth must be greater than or equal to one")
        sys.exit(1)

    # sanity check!
    if depth == 1:
        logger.info("Truncation depth of 1 may not be very useful")

    numprocs = op.numprocs
    if numprocs is not None and numprocs < 1:
        #add automatic parallelisation
        if numprocs==0:
            numprocs = mp.cpu_count()
        else:
            sys.stderr.write("Error: can't have a negative worker process count")
            sys.exit(1)

    if op.gpprior is None or op.gpprior == 'uniform':
        gpprior = None
    else:
        try:
            gpprior = parse_gp_hyperparam_priors(op.gpprior)
        except ValueError(s):
            sys.stderr.write("Error: "+s)
            sys.exit(1)

    # figure out where our output is going
    if op.csvoutput is None:
        csvoutput = None
    else:
        if op.csvoutput == '-':
            fd = sys.stdout
        else:
            fd = open(op.csvoutput,'w')
        csvoutput = csv.writer(fd)

    if op.hdf5output:
        hdf5output = h5.File(op.hdf5output,'w')
    else:
        hdf5output = None

    if hdf5output is None and csvoutput is None:
        logger.warning("No output will be saved, "
                       "this is only useful for debugging and benchmarking.")

    # load the data from disk
    inp = csi.loadData(fname[0])

    # check whether the second level is sorted (currently check whether all
    # levels are sorted, need to fix!)
    assert (inp.columns.is_monotonic_increasing)
    # not sure whether I can do anything similar for the rows
    
    #normalise the data
    if op.normalise == 'standardise':
        inp[:][:] = sp.stats.mstats.zscore(inp,axis=1,ddof=1)
    elif op.normalise == 'center':
        inp[:][:] = inp[:][:] - np.mean(inp[:][:],axis=1)[:,None]

    if op.verbose:
        logger.info("Genes: %s",
                    ", ".join([repr(x) for x in inp.index]))
        logger.info("Treatments: %s",
                    ", ".join([repr(x) for x in inp.columns.levels[0]]))
        logger.info("Time: %s",
                    ", ".join([repr(x) for x in inp.columns.levels[1]]))
        if gpprior is None:
            logger.info("Hyperparameters: uniform")
        else:
            logger.info("Hyperparameters: Gamma({0},{1})".format(*gpprior))

    # figure out which genes/rows we're going to process
    genes = op.genes
    if genes is None:
        logger.debug("No genes specified, assuming all")
        genes = list(inp.index)
    else:
        missing = np.setdiff1d(genes, inp.index)
        if len(missing) > 0:
            sys.stderr.write("Error: The following genes were not found: {missing}\n".format(
                missing=', '.join(missing)))
            sys.exit(1)

    # TODO: how does the user specify the parental set?
    if op.tflist is not None:
        with open(op.tflist,'r') as fid:
            tfs = fid.readlines()
        #lose the new line symbol at the end
        tfs = [x.rstrip() for x in tfs]
        missing = np.setdiff1d(tfs, inp.index)
        if len(missing) > 0:
            sys.stderr.write("Error: The following transcription factors were not found: {missing}\n".format(
                missing=', '.join(missing)))
            sys.exit(1)
    else:
        #assume everything is TFs
        logger.debug("No transcription factor list specified, assuming all")
        tfs = list(inp.index) 

    cc = csi.Csi(inp, tfs)
    em = cc.getEm()

    if gpprior:
        em.set_priors(gpprior[0], gpprior[1])

    if hdf5output:
        cc.write_hdf5(hdf5output)
        hdf5output.flush()

    if op.weighttrunc:
        val = float(op.weighttrunc)
        if not (0 < val < 1):
            sys.stderr.write("Error: The weight truncation must be between zero and one\n")
            sys.exit(1)

        if val > 0.01:
            logger.warning("weight truncation should probably be less than 0.01")

        em.weighttrunc = val

    if op.initweights:
        if op.initweights == 'uniform':
            em.sampleinitweights = False
        elif op.initweights == 'weighted':
            em.sampleinitweights = True
        else:
            sys.stderr.write("Error: Unrecognised initial weight mode: {initweights}\n".join(
                initweights=op.initweights))
            sys.exit(1)

    for i,res in enumerate(csi.runCsiEm(em, genes, lambda gene: cc.allParents(gene,depth), numprocs)):
        if csvoutput:
            res.writeCsv(csvoutput)
        if hdf5output:
            res.write_hdf5(hdf5output, i)
            hdf5output.flush()

if __name__ == '__main__':
    main()
