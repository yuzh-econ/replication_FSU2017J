import numpy as np
import pandas as pd
from numpy.linalg import inv

"""
#####################
##  Original Code  ##
#####################
"""
def MedianAD(arr, axis=None):
    """
    Compute Median Absolute Deviation (MAD) of an array along the given axis.

    Parameters
    ----------
    arr : array_like
        Input array or object that can be converted to an array.
    axis : int or None, optional
        Axis along which to compute MAD. If None, compute MAD over the flattened array.

    Returns
    -------
    mad : ndarray
        Median Absolute Deviation of the input array.
    """
    isDF = type(arr) is pd.DataFrame
    cols = None
    index = None

    if isDF:
        cols = list(arr)
        index = list(arr.index)
        arr = np.array((arr))

    # Compute median along the specified axis
    median = np.nanmedian(arr, axis=axis, )

    # Compute absolute deviations from the median
    absdev = np.abs(arr - median)

    # Compute median of absolute deviations along the specified axis
    mad = np.nanmedian(absdev, axis=axis)

    if isDF:
        if axis is None:
            return mad
        elif axis == 0:
            return pd.Series(mad, index=cols)
        else:
            return pd.Series(mad, index=index)

    else:
        return mad

def simulate_AR1(g, h, eta, T=10, x0=None, e=None):
    """
    Simulats time series from the VAR(1) model:
        x_t+1 = h x_t + eta e_t+1     for t=1,...,T-1
        y_t   = g x_t                 for t=1,...,T

    Parameters:
        g (ndarray): n_x by n_x matrix. Contemporaneous coefficient.
        h (ndarray): n_x by n_x matrix. Lagged coefficient.
        eta (ndarray): n_x by n_x matrix. Shock loading (coefficient).

        T (integer): Length of simulation. (Optional) Default = 250.
        x0 (ndarray): T by n_x matrix. Initial condition of x. (Optional) Default = 0.
        e (ndarray): T by n_e matrix. Random Shocks/ Innovation Process. (Optional) Defalut ~ iid N(0,I).

    Returns:
        X (ndarray): T by n_x matrix. Time series for variable X.
        Y (ndarray): T by n_y matrix. Time series for variable Y.
        e (ndarray): T by n_e matrix. Random Shocks used.
    """
    nx = h.shape[0]  # number of variables in X
    ny = g.shape[0]  # number of variables in Y
    ne = eta.shape[1]  # number of shocks

    # Initialize innovation process
    if e is None:
        e = np.random.randn(T, ne)
    elif e.shape != (T, ne):
        if e.shape[0] != T:
            T = e.shape[0]
        if e.shape[1] != ne:
            raise ValueError("Check dimension of errors.")

    # Initialize X and Y
    X = np.zeros([T, nx])
    Y = np.zeros([T, ny])

    # Initialize x0
    if x0 is None:
        x0 = np.zeros([nx, 1])
    elif x0.shape != ([nx, 1]):
        raise ValueError("Check dimensions of x0.")

    # Intialize Time Series
    X[0, :] = x0.T
    Y[0, :] = np.matmul(X[0, :], g.T)

    # Propagate Time Series
    for t in range(1, T):
        X[t, :] = np.matmul(X[t-1, :], h.T) + np.matmul(e[t, :], eta.T)
        Y[t, :] = np.matmul(X[t, :], g.T)

    return Y, X, e

def estimate_bias(hx, ETA, Tshort, Tlong, npr, Tmontecarlo):
    """
    Estimate the small sample bias of a structural VAR using Monte Carlo simulation.

    Parameters
    ----------
    hx : array-like, shape (n_vars, n_vars)
        Autoregressive matrices for the structural VAR.
    ETA : array-like, shape (n_vars, n_obs)
        Residuals of the structural VAR.
    Tshort : int
        Length of the short sample used to estimate the domestic block.
    Tlong : int
        Length of the long sample used to estimate the commodity price block.
    npr : int
        Number of variables in the commodity price block.
    Tmontecarlo : int
        Number of Monte Carlo simulations to perform.

    Returns
    -------
    sB : float
        Estimated small sample bias.
    stdB : float
        Standard deviation of the estimated small sample bias.
    Shx : array-like, shape (n_vars, n_vars)
        Average autoregressive matrices estimated over the Monte Carlo simulations.
    SETA : array-like, shape (n_vars, n_vars)
        Average covariance matrix of the residuals estimated over the Monte Carlo simulations.
    """
    nv = hx.shape[0] # number of variables in the SVAR
    gx = np.eye(nv) # identify matrix

    # True variance decomposition
    trueVD = variance_decomposition(gx, hx, ETA)[0]
    trueVS = trueVD[:npr, -(nv-npr):].sum(axis=0)
    # trueVS = np.sum(trueVD[0:npr, npr:nv])

    vsVS = np.empty([1,nv-npr])
    vsVS[:] = np.nan
    # vsVS = []
    Shx = np.zeros_like(hx)
    SETA = np.zeros_like(ETA)

    for i in range(1,Tmontecarlo+1):
        # Simulate artificial time series
        Y, _X, _e = simulate_AR1(gx, hx, ETA, 250)

        # Use the long sample to estimate the commodity price block
        # and the short sample to estimate the domestic block.
        Dlong = Y[-Tlong:]
        Dshort = Y[-Tshort:]

        # Estimate commodity price block
        p = Dlong[1:, 0:npr]
        p1 = np.roll(Dlong[:, 0:npr], 1, axis=0)[1:, 0:npr]
        T = p.shape[0]
        cons = np.ones((T, 1))
        X = np.hstack([p1, cons])
        b = inv(X.T @ X) @ X.T @ p
        A = b[:-1].T
        mu = p - X @ b
        Sigma_mu = np.cov(mu.T)

        # Estimate domestic block
        D = Dshort[1:, :]
        D1 = np.roll(Dshort, 1, axis=0)[1:, :]  # todo: shouldbeok
        T = D.shape[0]
        cons = np.ones((T, 1))
        X = np.hstack([D1, D[:,:npr], cons])
        Y = Dshort[1:, npr:nv]
        b_cons = inv(X.T @ X) @ X.T @ Y
        b_noc = b_cons[:-1]
        e = Y - X @ b_cons
        B = b_noc[:npr].T
        C = b_noc[npr:nv].T
        D = b_noc[-npr:].T
        Sigma_eps = np.cov(e.T)

        # Combine blocks to get autoregressive matrix for the structural VAR
        shx = np.block([[A, np.zeros((npr, nv-npr))],
                    [D @ A + B, C]])
        G = np.block([[np.eye(npr), np.zeros((npr, nv-npr))], [D, np.eye(nv-npr)]])
        Sigma = np.block([[Sigma_mu, np.zeros((npr, nv-npr))], [np.zeros((nv-npr, npr)), Sigma_eps]])

        # variance decomposition
        sETA = np.linalg.cholesky(G @ Sigma @ G.T)

        # To filter a possible small sample bias in the estimation of the price VAR, unpercentage the following two lines:
        # shx[:np, :np] = hxp
        # sETA[:np, :np] = ETAp

        V = variance_decomposition(gx, shx, sETA)[0]

        vsVS = np.vstack([vsVS, V[:npr,-(nv-npr):].sum(axis=0)])
        # vsVS = np.concatenate([vsVS, np.sum(V[:npr, npr:], axis=0)])

        Shx = (i-1)/i * Shx + shx / i
        SETA = (i-1)/i * SETA + sETA / i

    sVS = np.nanmean(vsVS, axis=0)  # Compute the mean of the variance share estimates from the Monte Carlo simulations

    if isinstance(sVS, float):
        sVS = np.array([sVS])

    if isinstance(trueVS, float):
        trueVS = np.array([trueVS])

    if sVS.shape != trueVS.shape:
        raise ValueError("sVS not computed correctly")  # Raise an error if the shape of the variance share estimates is incorrect

    sB = sVS - trueVS       # Compute the bias estimate by subtracting the true variance share from the mean variance share
    stdB = np.nanstd(sVS)   # Compute the standard deviation of the variance share estimates from the Monte Carlo simulations

    return sB, stdB, Shx, SETA  # Return the bias estimate, standard deviation, estimated VAR coefficients and variance-covariance matrix





"""
#####################
## Translated Code ##
#####################
"""


def variance_decomposition(gx, hx, ETA1, ETA2=None):
    """
    Computes the variance decomposition of y and x given matrices gx, hx, ETA1, and ETA2.
    If ETA2 is not provided, only the variance decomposition of y and x given ETA1 is computed.
    
    Parameters:
        gx (ndarray): n_y x n_x matrix.
        hx (ndarray): n_x x n_x matrix.
        ETA1 (ndarray): n_x x n_eta1 matrix.
        ETA2 (ndarray): n_x x n_eta2 matrix (optional).
    
    Returns:
        Vyr (ndarray): n_eta1 x n_y matrix of variance decompositions of y_t given ETA1.
        Vxr (ndarray): n_eta1 x n_x matrix of variance decompositions of x_t given ETA1.
        Vy (ndarray): n_y x n_y matrix of variance decompositions of y_t.
        Vx (ndarray): n_x x n_x matrix of variance decompositions of x_t.
    """
    
    Vy = []
    Vx = []
    
    # Compute variance decomposition given ETA1
    n1 = ETA1.shape[1]
    for j in range(n1):
        I1 = np.zeros((n1, n1))
        I1[j, j] = 1
        
        V1 = ETA1 @ I1 @ ETA1.T
        sigy, sigx = mom(gx, hx, V1)
        
        Vy.append(np.diag(sigy))
        Vx.append(np.diag(sigx))
        
    Vy = np.array(Vy)
    Vx = np.array(Vx)
    
    # Compute variance decomposition given ETA2, if provided
    if ETA2 is not None:
        n2 = ETA2.shape[1]
        for j in range(n2):
            I2 = np.zeros((n2, n2))
            I2[j, j] = 1
            
            V2 = ETA2 @ I2 @ ETA2.T
            sigy, _ = mom(gx, hx, V2)
            
            Vy = np.vstack((Vy, np.diag(sigy)))
    
    # Normalize the matrices
    Vyr = Vy / np.sum(Vy, axis=0)
    Vxr = Vx / np.sum(Vx, axis=0)
    
    return Vyr, Vxr, Vy, Vx


def mom(gx, hx, varshock, J=0, method=1):
    """
    Computes the unconditional variance-covariance matrix of x(t) with x(t+J),
    that is sigxJ=E[x(t)*x(t+J)'], and the unconditional variance covariaance matrix
    of y(t) with y(t+J), that is sigyJ=E[y(t)*y(t+J)'] where x(t) evolves as
    x(t+1) = hx x(t) + e(t+1) and y(t) evolves according to y(t) = gx x(t)
    where E[e(t)e(t)']=varshock. The parameter J can be any integer.
    method=1: use doubling algorithm, method!=1: use algebraic method.
    :param gx: array_like
        Coefficient matrix for the observation equation y(t) = gx x(t).
    :param hx: array_like
        Coefficient matrix for the state equation x(t+1) = hx x(t) + e(t+1).
    :param varshock: array_like
        Covariance matrix of the state disturbance e(t+1).
    :param J: int, optional
        Time horizon.
    :param method: int, optional
        Algorithm for computing the variance-covariance matrices.
    :return:
        sigyJ : ndarray
            The unconditional variance-covariance matrix of y(t) and y(t+J).
        sigxJ : ndarray
            The unconditional variance-covariance matrix of x(t) and x(t+J).
    """

    if J == 0:
        # If J is not provided, set it to 0.
        pass

    if method == 1:
        # Doubling algorithm
        hx_old = hx
        sig_old = varshock
        sigx_old = np.eye(hx.shape[0])
        diferenz = 0.1
        while diferenz > 1e-25:
            sigx = hx_old @ sigx_old @ hx_old.T + sig_old
            diferenz = np.max(np.abs(sigx - sigx_old))
            sig_old = hx_old @ sig_old @ hx_old.T + sig_old
            hx_old = hx_old @ hx_old
            sigx_old = sigx.copy()

    else:
        # Algebraic method
        # Get the variance of x
        # sigx = inv(I - hx⊗hx) vec(varshock)
        kronecker = np.kron(hx, hx)
        sigx = np.linalg.solve(np.eye(hx.shape[0]**2) - kronecker, varshock.ravel()).reshape(hx.shape)

    if J != 0:
        # Get E{x(t)*x(t+J)'}
        sigxJ = hx ** (-min(0, J)) @ sigx @ hx.T ** (max(0, J))
    else:
        sigxJ = sigx

    # Get E{y(t)*y(t+J)'}
    sigyJ = gx @ sigxJ @ gx.T.conj().T.real

    return sigyJ, sigxJ





if __name__=='__main__':
    import pickle
    with open('filename.pickle', 'rb') as handle:
        b = pickle.load(handle)
    _totalT_, _F_, _eta_ = b['test']
    k = estimate_bias(_F_, _eta_, _totalT_, 55, 4, 1000)
    print(k[0])
