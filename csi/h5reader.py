import numpy as np
import h5py as h5

import gp

class Model(object):
    def __init__(self, h5fd):
        self.fd = h5fd

        self.items = [s.decode('utf-8') for s in h5fd['items']]
        self.reps  = [Replicate(d) for d in h5fd['data'].values()]

    def get_res(self, res):
        return Result(self, self.fd["result/{0}".format(res+1)])

    def iter_res(self):
        for res in self.fd["result"].values():
            yield Result(self, res)

class Replicate(object):
    def __init__(self, h5data):
        self.name = h5data.attrs["replicate"]
        self.time = h5data.attrs["time"]
        self.data = h5data[:]

class Result(object):
    def __init__(self, mod, res):
        self.mod = mod
        self.res = res

        self.hypers = self.res.attrs['hyperparams']
        self.target = self.res.attrs['item'][0]

        self.ll      = self.res["loglik"][:]
        self.psets   = self.res['parents'][:]
        self.weights = self.res['weight'][:]

    def filter_by_weight(self, min_weight):
        keep = np.nonzero(self.weights >= min_weight)
        self.ll      = self.ll[keep]
        self.psets   = self.psets[keep]
        self.weights = self.weights[keep]

    def sort_psets(self):
        ii = np.argsort(-self.ll)

        self.ll      = self.ll[ii]
        self.psets   = self.psets[ii]
        self.weights = self.weights[ii]

class Predictor(object):
    def __init__(self, res, pset, datasets=None):
        self.target = res.target
        self.hypers = res.hypers
        self.pset   = list(pset)
        self.ix     = [self.target]+self.pset
        self.iy     = [self.target]

        if datasets is None:
            datasets = (r.data for r in res.mod.reps)

        X = [] # parent array
        Y = [] # target array
        for d in datasets:
            X.append(d[self.ix,:-1])
            Y.append(d[self.iy,1:])
        X = np.hstack(X).T
        Y = np.hstack(Y).T
        self.gp = gp.rbf(X,Y,self.hypers)

    def predict_dataset(self, data):
        return self.gp.predict(data[self.ix,:].T)

    def predict1(self, expr):
        return self.gp.predict(expr[None,self.ix])
