import numpy as np

import scipy.spatial.distance as dist
import scipy.linalg.lapack as lapack

_LOG_2_PI = 1.837877066409345483556

# my attempt at an optimised version of the RBF kernel and gaussian
# likelihood from the GPy library.  All I care about is getting at the
# log-likelihood and gradient of parameters quickly for a given theta

class rbf(object):
    def __init__(self, X, Y, theta):
        theta = theta ** [2,1,2]
        [sf2, l2, sn2] = theta

        # evaluate RBF kernel for our given X
        r = dist.pdist(X) / l2
        K = dist.squareform(sf2 * np.exp(-0.5 * r**2))
        np.fill_diagonal(K, sf2)

        # add in Gaussian noise (+ a bit for numerical stability)
        Ky = K.copy()
        np.fill_diagonal(Ky, sf2 + sn2 + 1e-8)

        # compute the Cholesky factorization of our covariance matrix
        LW, info = lapack.dpotrf(Ky, lower=True)
        assert info == 0

        # calculate lower half of inverse of K (assumes real symmetric positive definite)
        Wi, info = lapack.dpotri(LW, lower=True)
        assert info == 0

        # make symmetric by filling in the upper half
        Wi += np.tril(Wi,-1).T

        # and solve
        alpha, info = lapack.dpotrs(LW, Y, lower=True)
        assert info == 0

        # save these for later
        self.X = X
        self.Y = Y
        self.theta = theta
        self.r = r
        self.K = K
        self.Ky = Ky
        self.LW = LW
        self.Wi = Wi
        self.alpha = alpha

    def log_marginal(self):
        X     = self.X
        Y     = self.Y
        LW    = self.LW
        alpha = self.alpha

        # this lets us get at the log-determinant…
        W_logdet = 2.*np.sum(np.log(np.diag(LW)))
        #   and hence to the log of the marginal likelihood
        return 0.5 * (-Y.size * _LOG_2_PI - Y.shape[1] * W_logdet - np.sum(alpha * Y))

    def gradient_theta(self):
        X     = self.X
        Y     = self.Y
        r     = self.r
        K     = self.K
        Wi    = self.Wi
        LW    = self.LW
        alpha = self.alpha

        [sf2, l2, sn2] = self.theta

        dL_dK  = 0.5 * (np.dot(alpha,alpha.T) - Y.shape[1] * Wi)

        gradient = np.zeros(3)
        # fill in gradient of sn2
        gradient[2] = np.diag(dL_dK).sum() * np.sqrt(sn2)*2

        # multiply to save duplication of work
        dL_dK *= K
        # gradient of sf2
        gradient[0] = np.sum(dL_dK) / sf2 * np.sqrt(sf2)*2
        # gradient of l2
        gradient[1] = np.sum(dist.squareform(r)**2 * dL_dK) / l2

        return gradient

    def predict(self, X2):
        X     = self.X
        Wi    = self.Wi
        alpha = self.alpha

        [sf2, l2, sn2] = self.theta

        r = dist.cdist(X,X2) / l2
        Kx = sf2 * np.exp(-0.5 * r**2)

        WiKx = np.dot(Wi, Kx)

        mu = np.dot(Kx.T, alpha)
        var = (sf2 - np.sum(WiKx*Kx, 0))[:,None]

        return (mu,var)
