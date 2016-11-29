import numpy as np
import scipy as sp
import scipy.optimize as spo
import scipy.stats as sps
import pandas as pd
import h5py as h5

import multiprocessing as mp

import itertools as it

import gp

import logging
logger = logging.getLogger('CSI')

def parentalSets(items, item, depth, tfs=None):
    """Iterate over all "Parental Sets".

    A parental set is a list of regulators/transcription factors for a
    specific gene.  The parental sets for any given item is every
    subset of items up to a given depth that does not include the
    item.
    """

    # filter to TFs if needed
    if tfs is not None:
        # duplicate list to avoid modifying callers state
        items = list(tfs)

    # exclude the target if needed
    if item in items:
        # duplicate list to avoid modifying callers state
        items = list(items)
        items.remove(item)

    for i in range(0, depth+1):
        # iterate over every subset of size i
        for subset in it.combinations(items, i):
            yield (list(subset),item)

def logexp_optimise(fn, x, **kwargs):
    """Our 'LogExp' transform is taken from GPy and ensures that
    parameters are always greater than zero.
    """
    def transform(x):
        theta = logexp_to_natural(x)
        y,grad = fn(theta)
        # get gradients back to our space
        grad *= logexp_gradientfactor(theta)
        return (y,grad)
    res = spo.minimize(transform, natural_to_logexp(x), jac=True,
                       **kwargs)
    res.x = logexp_to_natural(res.x)
    return res

def logexp_to_natural(x):
    return np.logaddexp(0., x)

def natural_to_logexp(f):
    # wrap in np.where to stop exp() of large values overflowing and
    # hence log giving back non-sensible values
    return np.where(f > 36., f, np.log(np.expm1(f)))

def logexp_gradientfactor(f):
    # no point wrapping in np.where as only small values will overflow
    # and they'll do that anyway
    return -np.expm1(-f)

class CsiError(Exception):
    pass

class CsiEmFailed(CsiError):
    def __init__(self, res):
        super(CsiEmFailed, self).__init__("Failed to optimise parameters ({msg})".format(
            msg=repr(res.message)))
        self.res = res

class CsiResult(object):
    def write_hdf5(self, file, num):
        raise NotImplementedError

class EmRes(CsiResult):
    def __init__(self, em):
        self.em      = em
        self.hypers  = np.copy(em.hypers)
        self.weights = np.copy(em.weights)
        self.pset    = em.pset
        self.ll      = np.copy(em.logliks())

    def writeCsv(self, csvout):
        for i in np.argsort(-self.ll):
            pset    = self.pset[i]
            weights = self.weights[i]
            csvout.writerow(
                [pset[1],
                 ':'.join(pset[0]),
                 weights]+
                self.hypers.tolist())

    def getMarginalWeights(self):
        # add all genes to our (directed) graph
        for pset,weight in zip(self.pset,self.weights):
            target = pset[1]
            for regulator in pset[0]:
                yield (regulator,target,weight)

    def enum_predictions(self):
        reps = self.em.csi.get_replicates()
        ilocs = [np.array(r.iloc[:-1]) for r in reps]

        for pset in self.pset:
            mu,var = self.em._predict_pset(pset, self.hypers)

            yield dict(
                mu  = [mu [l-i,0].tolist() for i,l in enumerate(ilocs)],
                var = [var[l-i,0].tolist() for i,l in enumerate(ilocs)]
            )

    def write_hdf5(self, file, num):
        itemmap = self.em.csi._itemmap
        ptype = h5.special_dtype(enum=('i', itemmap))

        pset = np.array([np.array([itemmap[i] for i in pi],dtype=np.int32)
                         for pi,gi in self.pset],
                        dtype=h5.special_dtype(vlen=ptype))

        grp = file.create_group("result/{n}".format(n=num+1))

        ga = grp.attrs
        ga['hyperparams'] = self.hypers
        ga['restype']     = 'EM'
        ga['item']        = np.array([itemmap[self.pset[0][1]]],
                                     dtype=ptype)

        grp.create_dataset('loglik',data=self.ll)
        grp.create_dataset('weight',data=self.weights)
        grp.create_dataset('parents',data=pset)

class CsiEmWorker(object):
    def __init__(self, X, Y, pset):
        self.X = X
        self.Y = Y
        self.pset = pset

    def eval_gp(self, i, theta):
        "Evaluate a GP for parental set @i using a RBF"

        # extract the items for our given parental set
        pi,gi = self.pset[i]

        # pull out data for these items
        X = self.X.loc[:,[gi]+pi].values
        Y = self.Y.loc[:,[gi]].values

        # evaluate the GP
        return gp.rbf(X, Y, theta)

    def loglik_grad(self, args):
        "Calculate marginal likelihood and gradient of given parameters"
        i,weight,theta = args

        m = self.eval_gp(i, theta)

        return (weight * m.log_marginal(), weight * m.gradient_theta())

    def loglik(self, args):
        "Calculate the marginal likelihood for our parental set"
        i,theta = args

        return self.eval_gp(i, theta).log_marginal()

def _pool_init(X, Y, pset):
    "Initialise (global) variables shared by all Pool workers"
    global _pool
    _pool = CsiEmWorker(X, Y, pset)

def _pool_loglik_grad(args):  return _pool.loglik_grad(args)
def _pool_loglik(args):       return _pool.loglik(args)

class CsiEm(object):
    "Find MAP estimate via Expectation-Maximisation."
    def __init__(self, csi):
        self.csi = csi
        self.X = csi.X
        self.Y = csi.Y
        self.pset = []
        self.weighttrunc = 1e-5
        self.sampleinitweights = True
        self._prior_shape = None
        self._prior_scale = None

        self.pool = None

    @property
    def hypers(self):
        return self._hypers

    @hypers.setter
    def hypers(self, value):
        self._hypers   = value
        self._updatell = True

    def set_priors(self, shape, scale):
        assert (shape is None) == (scale is None)
        if shape is None:
            self._prior_shape = None
            self._prior_scale = None
            return
        shape = np.array(shape)
        scale = np.array(scale)
        assert shape.shape in [(),(3,)]
        assert scale.shape in [(),(3,)]
        self._prior_shape = shape
        self._prior_scale = scale

    def setup(self, pset, poolsize=None):
        "Configure model for EM using the specified parent set."

        if self.sampleinitweights:
            # want to start with mostly singleton parental sets, so start
            # by calculating parental set sizes
            pl = np.array([len(a) for a,b in pset])
            # down-weight higher order parental sets
            w = sps.gamma.rvs(np.where(pl > 1, 1e-2, 0.5))
            # normalise weights to sum to one
            w /= np.sum(w)
        else:
            # all weights are equal
            w = np.ones(len(pset)) / len(pset)

        self.pset = pset
        self.hypers = sp.exp(sp.randn(3))
        self.weights = w
        self.tol = 0.1

        if poolsize is None:
            self.pool = None
            self.worker = CsiEmWorker(self.X, self.Y, self.pset)
        else:
            if self.pool is not None:
                self.pool.terminate()
            self.pool = mp.Pool(poolsize, _pool_init,
                                (self.X, self.Y, self.pset))

    def _predict_pset(self, pset, theta):
        i,j = pset

        X = self.X.loc[:,[j]+i].values
        Y = self.Y.loc[:,[j]].values

        return gp.rbf(X,Y,theta).predict(X)

    def _optfn(self, x):
        """Return negated-loglik and gradient in a form suitable for use with
        SciPy's numeric optimisation."""

        logger.debug("     optfn(theta=%s)", str(x))

        wmx = max(self.weights) * self.weighttrunc

        ip = []
        for i,w in enumerate(self.weights):
            if w < wmx:
                continue
            ip.append((i,w,x))

        if self.pool is None:
            itr = map(self.worker.loglik_grad, ip)
        else:
            itr = self.pool.imap_unordered(_pool_loglik_grad, ip, 10)

        if self._prior_shape is None:
            ll = 0.
            grad = np.zeros(len(x))
        else:
            ll   = sum(sp.special.xlogy(self._prior_shape-1,x)-(x/self._prior_scale))
            grad = (self._prior_shape - 1)/x - 1/self._prior_scale

        for l,g in itr:
            ll   += l
            grad += g

        logger.debug("       optfn=%g", ll)

        return -ll, -grad

    def optimiseHypers(self):
        """Re-optimise hyper-parameters of the given model, throwing exception
        on failure.
        """
        res = logexp_optimise(self._optfn, self.hypers, tol=self.tol)
        if not res.success:
            raise CsiEmFailed(res)
        self.hypers = res.x
        return res

    def logliks(self):
        "Return the log-likelihood of each GP given the current hyperparameters"
        if not self._updatell:
            return self._ll

        ip = zip(range(len(self.pset)),it.repeat(self.hypers))
        if self.pool is None:
            itr = list(map(self.worker.loglik, ip))
        else:
            itr = self.pool.map(_pool_loglik, ip, 10)

        self._ll = np.array(itr)
        self._updatell = False

        return self._ll

    def reweight(self):
        "Recalculate the weights and return the KL divergence."
        ll = self.logliks()
        # calculate the new weights
        w = np.exp(ll - max(ll))
        w /= np.sum(w)
        # update weights, saving old so we can calculate the KL divergence
        w0 = self.weights
        self.weights = w
        # calculate the KL divergence
        kl = sp.special.xlogy(w,w / w0).sum()

        # update optimiser tolerance
        self.tol = min(self.tol,max(kl/3,1e-4))
        return kl

    def getResults(self):
        return EmRes(self)

class CsiRep(object):
    def __init__(self,name,ilocs,time):
        self.name = name
        self.iloc = ilocs
        self.time = time

class Csi(object):
    "Perform CSI analysis on the provided data"
    def __init__(self, data, tfs=None):
        self.data = data
        self._items = data.index.tolist()
        self._itemmap = dict((b,a) for (a,b) in enumerate(self._items))
        self.tfs = tfs

        n = [a for a,b in data.columns.values]
        ix = np.flatnonzero(np.array([a==b for a,b in zip(n, n[1:])]))
        self.X = data.iloc[:,ix].T
        self.Y = data.iloc[:,ix+1].T

    def get_replicates(self):
        ""
        def c1(tup):
            a,(b,c) = tup
            return b
        def mk(a,l):
            l  = list(l)
            ix = [i for i,_ in l]
            v  = [v for _,(_,v) in l]
            return CsiRep(a,ix,v)

        return [mk(a,l) for a,l in it.groupby(enumerate(self.data.columns),c1)]

    def get_items(self):
        return self._items

    def item_to_num(self, item):
        return self._itemmap[item]

    def allParents(self, item, depth):
        "Utility function to calculate the parental set of the given item"
        return list(parentalSets(self._items, item, depth, self.tfs))

    def getEm(self):
        "For getting at a MAP estimate via expectation-maximisation."
        return CsiEm(self)

    def write_hdf5(self, file):
        file.create_dataset('items', data=np.string_(self._items))

        data = file.create_group('data')
        for i,r in enumerate(self.get_replicates()):
            d = self.data.iloc[:,r.iloc]
            df = data.create_dataset(str(i+1),data=d.values)
            df.attrs["replicate"] = r.name
            df.attrs["time"] = r.time

def loadData(path):
    with open(path) as fd:
        # read in data and make sure headers are correct
        inp = pd.read_csv(fd,dtype=str,index_col=0,header=[0,1])
        inp.columns = pd.MultiIndex.from_tuples([(a,float(b)) for a,b in inp.columns],
                                                names=inp.columns.names)
        # convert to floating point values
        return inp.astype(float)

def runCsiEm(em, genes, fnpset, poolsize=None):
    for gene in genes:
        logger.info("Processing: %s", repr(gene))

        em.setup(fnpset(gene), poolsize)

        for ittr in range(1, 20):
            logger.debug("%2i: optimising hyperparameters (tol=%g)",
                         ittr, em.tol)
            em.optimiseHypers()
            logger.debug("%2i: recalculating weights",
                         ittr)
            kl = em.reweight()
            logger.info("%2i: kl=%10.4g, hypers=%s",
                         ittr,kl,str(em.hypers))
            if kl < 1e-5:
                break

        logger.debug("finished after %i iterations", ittr)

        yield em.getResults()
