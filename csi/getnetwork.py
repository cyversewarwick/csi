import json
import sys

import pandas as pd
import h5py as h5
import h5reader as pp

import optparse
import logging

def update_marginal(arr, mod):
    for res in mod.iter_res():
        for w,ps in zip(res.weights,res.psets):
            arr[res.target,ps] += w

def update_map(arr, mod):
    for res in mod.iter_res():
        res.sort_psets()
        if len(res.weights) > 0:
            arr[res.target,res.psets[0]] += res.weights[0]

def cmdparser(args):
    # create parser objects
    op = optparse.OptionParser()
    out = optparse.OptionGroup(op, "Output Formats")

    op.set_usage("usage: %prog [options] FILE1.hdf5 [FILE2.hdf5]")
    op.set_defaults(
        verbose=False,
        mode=update_marginal)

    op.add_option('--map',action='store_const', const=update_map, dest='mode',
                  help="Generate MAP network")

    op.add_option('--marginal',action='store_const', const=update_marginal, dest='mode',
                  help="Generate marginal network")

    # define general parameters
    op.add_option('-v','--verbose',dest='verbose',action='count',
                  help="Increase verbosity (specify twice for more detail)")

    # parse our command line arguments
    return op.parse_args(args)

def main(args = None):
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

    if len(fname) < 1:
        sys.stderr.write("Error: I need at least one file to process,\n"
                         "  rerun with --help to get a summary of accepted options.\n")
        sys.exit(1)

    updatefn = op.mode

    try:
        with h5.File(fname[0],"r") as fd:
            mod = pp.Model(fd)
            df = pd.DataFrame(0,dtype=float,index=mod.items,columns=mod.items)
            arr = df.values
            updatefn(arr,mod)
    except IOError as err:
        sys.stderr.write("Error opening file {path}: {err}\n".format(path=repr(fname[0]), err=str(err)))
        sys.exit(1)

    for path in fname[1:]:
        with h5.File(path,"r") as fd:
            updatefn(arr,pp.Model(fd))

    df[:] = arr
    df.to_csv(sys.stdout)

if __name__ == '__main__':
    main()
