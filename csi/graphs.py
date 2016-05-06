import math

# see http://www.johndcook.com/standard_deviation.html
# and http://www.johndcook.com/skewness_kurtosis.html
class RunningStats(object):
    def __init__(self):
        self._n = 0
        self._m = 0.0
        self._s = 0.0

    def push(self, x):
        self._n += 1
        if self._n == 0:
            self._m = float(x)
            self._s = 0
        else:
            dx = x - self._m
            self._m += dx / self._n;
            self._s += dx * (float(x) - self._m);

    def __str__(self):
        return ("RunningStats(%i,%g,%g)" %
                (self.numel(),
                 self.mean(),
                 self.std()))

    def numel(self):
        return self._n

    def mean(self):
        return self._m

    def var(self):
        return self._s / (self._n-1)

    def std(self):
        return math.sqrt(self.var())


class CsiGraph(object):
    def __init__(self):
        self._g = {}

    def push(self, regulator, target, weight):
        key = (regulator,target)
        if key in self._g:
            v = self._g[key]
        else:
            v = RunningStats()
            self._g[key] = v
        v.push(weight)

    def __iter__(self):
        for k,v in self._g.items():
            yield (k[0],k[1],v.mean())
