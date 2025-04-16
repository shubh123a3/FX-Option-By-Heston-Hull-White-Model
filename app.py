import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy.special as sp
import scipy.integrate as integrate
import scipy.optimize as optimize
import enum
import numpy as np
import streamlit as st
class OptionType(enum.Enum):

    CALL = 1.0
    PUT = -1.0

 #Black-Scholes Call option price
def BS_Call_Put_Option_Price(CP,S_0,K,sigma,tau,r):
    if K is list:
        K = np.array(K).reshape([len(K),1])
    d1    = (np.log(S_0 / K) + (r + 0.5 * np.power(sigma,2.0))
    * tau) / (sigma * np.sqrt(tau))
    d2    = d1 - sigma * np.sqrt(tau)
    if CP == OptionType.CALL:
        value = stats.norm.cdf(d1) * S_0 - stats.norm.cdf(d2) * K * np.exp(-r * tau)
    elif CP == OptionType.PUT:
        value = stats.norm.cdf(-d2) * K * np.exp(-r * tau) - stats.norm.cdf(-d1)*S_0
    return value

# Implied volatility method
def ImpliedVolatilityBlack76(CP,marketPrice,K,T,S_0):
    # To determine initial volatility we interpolate define a grid for sigma
    # and interpolate on the inverse
    sigmaGrid = np.linspace(0.0,5.0,5000)
    optPriceGrid = BS_Call_Put_Option_Price(CP,S_0,K,sigmaGrid,T,0.0)
    sigmaInitial = np.interp(marketPrice,optPriceGrid,sigmaGrid)
    print("Strike = {0}".format(K))
    print("Initial volatility = {0}".format(sigmaInitial))

    # Use determined input for the local-search (final tuning)
    func = lambda sigma: np.power(BS_Call_Put_Option_Price(CP, S_0, K, sigma, T, 0.0) - marketPrice, 1.0)
    impliedVol = optimize.newton(func, sigmaInitial, tol=1e-15)
    print("Final volatility = {0}".format(impliedVol))
    if impliedVol > 2.0:
        impliedVol = 0.0
    return impliedVol

def CallPutOptionPriceCOSMthd_StochIR(cf,CP,S0,tau,K,N,L,P0T):
    # cf   - characteristic function as a functon, in the book denoted as \varphi
    # CP   - C for call and P for put
    # S0   - Initial stock price
    # tau  - time to maturity
    # K    - list of strikes
    # N    - Number of expansion terms
    # L    - size of truncation domain (typ.:L=8 or L=10)
    # P0T  - zero-coupon bond for maturity T.
    if K is not np.array:
        K = np.array(K).reshape([len(K), 1])

        # assigning i=sqrt(-1)
    i = 1j
    x0 = np.log(S0 / K)

    # truncation domain
    a = 0.0 - L * np.sqrt(tau)
    b = 0.0 + L * np.sqrt(tau)

    # sumation from k = 0 to k=N-1
    k = np.linspace(0, N - 1, N).reshape([N, 1])
    u = k * np.pi / (b - a)
    H_k = CallPutCoefficients(OptionType.PUT, a, b, k)
    mat = np.exp(i * np.outer((x0 - a), u))
    temp = cf(u) * H_k
    temp[0] = 0.5 * temp[0]
    value = K * np.real(mat.dot(temp))

    # we use call-put parity for call options
    if CP == OptionType.CALL:
        value = value + S0 - K * P0T

    return value


# Determine coefficients for Put Prices
def CallPutCoefficients(CP, a, b, k):
    if CP == OptionType.CALL:
        c = 0.0
        d = b
        coef = Chi_Psi(a, b, c, d, k)
        Chi_k = coef["chi"]
        Psi_k = coef["psi"]
        if a < b and b < 0.0:
            H_k = np.zeros([len(k), 1])
        else:
            H_k = 2.0 / (b - a) * (Chi_k - Psi_k)
    elif CP == OptionType.PUT:
        c = a
        d = 0.0
        coef = Chi_Psi(a, b, c, d, k)
        Chi_k = coef["chi"]
        Psi_k = coef["psi"]
        H_k = 2.0 / (b - a) * (- Chi_k + Psi_k)

    return H_k


def Chi_Psi(a, b, c, d, k):
    psi = np.sin(k * np.pi * (d - a) / (b - a)) - np.sin(k * np.pi * (c - a) / (b - a))
    psi[1:] = psi[1:] * (b - a) / (k[1:] * np.pi)
    psi[0] = d - c

    chi = 1.0 / (1.0 + np.power((k * np.pi / (b - a)), 2.0))
    expr1 = np.cos(k * np.pi * (d - a) / (b - a)) * np.exp(d) - np.cos(k * np.pi
                                                                       * (c - a) / (b - a)) * np.exp(c)
    expr2 = k * np.pi / (b - a) * np.sin(k * np.pi *
                                         (d - a) / (b - a)) - k * np.pi / (b - a) * np.sin(k
                                                                                           * np.pi * (c - a) / (
                                                                                                       b - a)) * np.exp(
        c)
    chi = chi * (expr1 + expr2)

    value = {"chi": chi, "psi": psi}
    return value

def EUOptionPriceFromMCPathsGeneralizedFXFrwd(CP,S,K):
    # S is a vector of Monte Carlo samples at T
    result = np.zeros([len(K),1])
    if CP == OptionType.CALL:
        for (idx,k) in enumerate(K):
            result[idx] = np.mean(np.maximum(S-k,0.0))
    elif CP == OptionType.PUT:
        for (idx,k) in enumerate(K):
            result[idx] = np.mean(np.maximum(k-S,0.0))
    return result


def GeneratePathsHHWFXHWEuler(NoOfPaths, NoOfSteps, T, frwdFX, v0, vbar, kappa, gamma, lambdd, lambdf, etad, etaf,
                              rhoxv, rhoxrd, rhoxrf, rhovrd, rhovrf, rhordrf):
    Wx = np.zeros([NoOfPaths, NoOfSteps + 1])
    Wv = np.zeros([NoOfPaths, NoOfSteps + 1])
    Wrd = np.zeros([NoOfPaths, NoOfSteps + 1])
    Wrf = np.zeros([NoOfPaths, NoOfSteps + 1])

    V = np.zeros([NoOfPaths, NoOfSteps + 1])
    FX = np.zeros([NoOfPaths, NoOfSteps + 1])
    V[:, 0] = v0
    FX[:, 0] = frwdFX

    dt = T / float(NoOfSteps)
    Bd = lambda t, T: 1.0 / lambdd * (np.exp(-lambdd * (T - t)) - 1.0)
    Bf = lambda t, T: 1.0 / lambdf * (np.exp(-lambdf * (T - t)) - 1.0)

    cov = np.array([[1.0, rhoxv, rhoxrd, rhoxrf], [rhoxv, 1.0, rhovrd, rhovrf], \
                    [rhoxrd, rhovrd, 1.0, rhordrf], [rhoxrf, rhovrf, rhordrf, 1.0]])

    time = np.zeros([NoOfSteps + 1])

    for i in range(0, NoOfSteps):
        Z = np.random.multivariate_normal([.0, .0, .0, .0], cov, NoOfPaths)
        if NoOfPaths > 1:
            Z[:, 0] = (Z[:, 0] - np.mean(Z[:, 0])) / np.std(Z[:, 0])
            Z[:, 1] = (Z[:, 1] - np.mean(Z[:, 1])) / np.std(Z[:, 1])
            Z[:, 2] = (Z[:, 2] - np.mean(Z[:, 2])) / np.std(Z[:, 2])
            Z[:, 3] = (Z[:, 3] - np.mean(Z[:, 3])) / np.std(Z[:, 3])

        Wx[:, i + 1] = Wx[:, i] + np.power(dt, 0.5) * Z[:, 0]
        Wv[:, i + 1] = Wv[:, i] + np.power(dt, 0.5) * Z[:, 1]
        Wrd[:, i + 1] = Wrd[:, i] + np.power(dt, 0.5) * Z[:, 2]
        Wrf[:, i + 1] = Wrf[:, i] + np.power(dt, 0.5) * Z[:, 3]

        # Variance process- Euler discretization
        V[:, i + 1] = V[:, i] + kappa * (vbar - V[:, i]) * dt \
                      + gamma * rhovrd * etad * Bd(time[i], T) * np.sqrt(V[:, i]) * dt \
                      + gamma * np.sqrt(V[:, i]) * (Wv[:, i + 1] - Wv[:, i])
        V[:, i + 1] = np.maximum(V[:, i + 1], 0.0)

        # FX process under the forward measure
        FX[:, i + 1] = FX[:, i] * (1.0 + np.sqrt(V[:, i]) * (Wx[:, i + 1] - Wx[:, i]) \
                                   - etad * Bd(time[i], T) * (Wrd[:, i + 1] - Wrd[:, i]) \
                                   + etaf * Bf(time[i], T) * (Wrf[:, i + 1] - Wrf[:, i]))
        time[i + 1] = time[i] + dt

    paths = {"time": time, "FX": FX}
    return paths


# Exact expectation E(sqrt(V(t)))
def meanSqrtV_3(kappa, v0, vbar, gamma):
    delta = 4.0 * kappa * vbar / gamma / gamma
    c = lambda t: 1.0 / (4.0 * kappa) * gamma * gamma * (1.0 - np.exp(-kappa * (t)))
    kappaBar = lambda t: 4.0 * kappa * v0 * np.exp(-kappa * t) / (gamma * gamma * (1.0 - np.exp(-kappa * t)))
    temp1 = lambda t: np.sqrt(2.0 * c(t)) * sp.gamma((1.0 + delta) / 2.0) / sp.gamma(delta / 2.0) * sp.hyp1f1(-0.5,
                                                                                                              delta / 2.0,
                                                                                                              -kappaBar(
                                                                                                                  t) / 2.0)
    return temp1


def C_H1HW_FX(u, tau, kappa, gamma, rhoxv):
    i = 1j

    D1 = np.sqrt(np.power(kappa - gamma * rhoxv * i * u, 2.0) + (u * u + i * u) * gamma * gamma)
    g = (kappa - gamma * rhoxv * i * u - D1) / (kappa - gamma * rhoxv * i * u + D1)
    C = (1.0 - np.exp(-D1 * tau)) / (gamma * gamma * (1.0 - g * np.exp(-D1 * tau))) \
        * (kappa - gamma * rhoxv * i * u - D1)
    return C

def GenerateStrikes(frwd,Ti):
    c_n = np.array([-1.5, -1.0, -0.5,0.0, 0.5, 1.0, 1.5])
    return frwd * np.exp(0.1 * c_n * np.sqrt(Ti))
def ChFH1HW_FX(u, tau, gamma, Rxv, Rxrd, Rxrf, Rrdrf, Rvrd, Rvrf, lambdd, etad, lambdf, etaf, kappa, vBar, v0):
    i = 1j
    C = lambda u, tau: C_H1HW_FX(u, tau, kappa, gamma, Rxv)
    Bd = lambda t, T: 1.0 / lambdd * (np.exp(-lambdd * (T - t)) - 1.0)
    Bf = lambda t, T: 1.0 / lambdf * (np.exp(-lambdf * (T - t)) - 1.0)
    G = meanSqrtV_3(kappa, v0, vBar, gamma)

    zeta = lambda t: (Rxrd * etad * Bd(t, tau) - Rxrf * etaf * Bf(t, tau)) * G(t) + \
                     Rrdrf * etad * etaf * Bd(t, tau) * Bf(t, tau) - 0.5 * (
                                 etad ** 2.0 * Bd(t, tau) ** 2.0 + etaf ** 2.0 * Bf(t, tau) ** 2.0)

    # Integration in the function A(u,tau)
    int1 = np.zeros([len(u), 1], dtype=np.complex_)
    N = 500
    z = np.linspace(0.0 + 1e-10, tau - 1e-10, N)

    temp1 = lambda z1: kappa * vBar + Rvrd * gamma * etad * G(tau - z1) * Bd(tau - z1, tau)
    temp2 = lambda z1, u: -Rvrd * gamma * etad * G(tau - z1) * Bd(tau - z1, tau) * i * u
    temp3 = lambda z1, u: Rvrf * gamma * etaf * G(tau - z1) * Bf(tau - z1, tau) * i * u
    f = lambda z1, u: (temp1(z1) + temp2(z1, u) + temp3(z1, u)) * C(u, z1)

    value1 = np.trapz(np.real(f(z, u)), z).reshape(u.size, 1)
    value2 = np.trapz(np.imag(f(z, u)), z).reshape(u.size, 1)
    int1 = (value1 + value2 * i)
    int2 = (u ** 2.0 + i * u) * np.trapz(zeta(tau - z), z)
    A = int1 + int2

    cf = np.exp(A + v0 * C(u, tau))
    return cf


def mainCalculation():
    st.title("FX Option Pricing with Hybrid Heston-Hull-white model")
    st.write("This example shows how to price FX options using the Heston-Hull-white model with the COS method.")
    st.write("The model is calibrated to the market data.")

    option_type = {"Call": OptionType.CALL, "Put": OptionType.PUT}
    slected_option = st.sidebar.selectbox("Select the option type", option_type.keys())
    CP = option_type[slected_option]
    T =st.sidebar.slider("Select the maturity", min_value=0.1, max_value=20.0, value=5.0, step=0.1)

    NoOfPaths = st.sidebar.slider("Select the number of paths", min_value=1000, max_value=100000, value=1000, step=1000)
    NoOfSteps = (int)(T * 50)

    # Settings for the COS method
    N = 500
    L = 8

    # Market Settings
    st.sidebar.write("Market settings")
    P0Td = lambda t: np.exp(-0.02 * t)
    P0Tf = lambda t: np.exp(-0.05 * t)
    y0 =st.sidebar.number_input("Select the initial forward FX rate y0", value=1.35, step=0.01)
    frwdFX = y0 * P0Tf(T) / P0Td(T)
    kappa = st.sidebar.number_input("Select the mean reversion speed kappa", value=0.5, step=0.01)
    gamma = st.sidebar.number_input("Select the volatility of variance gamma", value=0.3, step=0.01)
    vbar = st.sidebar.number_input("Select the long term variance vbar", value=0.1, step=0.01)
    v0 = st.sidebar.number_input("Select the initial variance v0", value=0.1, step=0.01)

    # HW model settings
    st.sidebar.write("Hull-white model settings")
    lambdd = st.sidebar.number_input("Select the mean reversion speed lambda_d domestic " , value=0.01, step=0.01)
    lambdf = st.sidebar.number_input("Select the mean reversion speed lambda_f forign process", value=0.05, step=0.01)
    etad = st.sidebar.number_input("Select the volatility of the domestic process", value=0.007, step=0.01)
    etaf = st.sidebar.number_input("Select the volatility of the forign process", value=0.012, step=0.01)

    # correlations
    st.sidebar.write("Correlations")
    Rxv = -0.4
    Rxrd = -0.15
    Rxrf = -0.15
    Rvrd = 0.3
    Rvrf = 0.3
    Rrdrf = 0.25
    st.write(""" Correlation domestic and forign process
    Rxv = -0.4
    Rxrd = -0.15
    Rxrf = -0.15
    Rvrd = 0.3
    Rvrf = 0.3
    Rrdrf = 0.25
    """)

    # Strikes
    K = GenerateStrikes(frwdFX, T)
    K = np.array(K).reshape([len(K), 1])

    # number of repeated simulations for different Monte Carlo seed
    SeedV = range(0, 20)
    optMCM = np.zeros([len(SeedV), len(K)])
    for (idx, seed) in enumerate(SeedV):
        print('Seed number = {0} out of {1}'.format(idx, len(SeedV)))
        np.random.seed(seed)
        paths = GeneratePathsHHWFXHWEuler(NoOfPaths, NoOfSteps, T, frwdFX, v0, vbar, kappa, gamma, lambdd, lambdf, etad,
                                          etaf, Rxv, Rxrd, Rxrf, Rvrd, Rvrf, Rrdrf)
        frwdfxT = paths["FX"]
        optMC = P0Td(T) * EUOptionPriceFromMCPathsGeneralizedFXFrwd(CP, frwdfxT[:, -1], K)
        optMCM[idx, :] = np.squeeze(optMC)

    # Average of the runs + standard deviation
    optionMC_E = np.zeros([len(K)])
    optionMC_StDev = np.zeros([len(K)])
    for (idx, k) in enumerate(K):
        optionMC_E[idx] = np.mean(optMCM[:, idx])
        optionMC_StDev[idx] = np.std(optMCM[:, idx])

    # Value from the COS method
    cf = lambda u: ChFH1HW_FX(u, T, gamma, Rxv, Rxrd, Rxrf, Rrdrf, Rvrd, Rvrf, lambdd, etad, lambdf, etaf, kappa, vbar,
                              v0)
    valCOS_H1HW = P0Td(T) * CallPutOptionPriceCOSMthd_StochIR(cf, CP, frwdFX, T, K, N, L, 1.0)

    # checking martingale property
    EyT = P0Td(T) / P0Tf(T) * EUOptionPriceFromMCPathsGeneralizedFXFrwd(CP, frwdfxT[:, -1], [0.0])
    print("Martingale check: P_d(T)/P_f(T)*E[FX(T)] ={0:.4f} and y0 ={1}".format(EyT[0][0], y0))
    print("Maturity chosen to T={0}".format(T))
    for (idx, k) in enumerate(K):
        st.success(
            "Option price for strike K={0:.4f} is equal to: COS method = {1:.4f} and MC = {2:.4f} with stdDev = {3:.4f}".format(
                k[0], valCOS_H1HW[idx][0], optionMC_E[idx], optionMC_StDev[idx]))

    fig1=plt.figure(1)
    plt.plot(K, optionMC_E)
    plt.plot(K, valCOS_H1HW, '--r')
    plt.grid()
    plt.legend(['Monte Carlo Option Price', 'COS method'])
    plt.title("Fx Option prices")
    plt.xlabel("Strike")
    plt.ylabel("Option price")
    st.pyplot(fig1)

    # implied volatilities
    IVCos = np.zeros([len(K), 1])
    IVMC = np.zeros([len(K), 1])
    for (idx, k) in enumerate(K):
        priceCOS = valCOS_H1HW[idx] / P0Td(T)
        IVCos[idx] = ImpliedVolatilityBlack76(CP, priceCOS, k, T, frwdFX) * 100.0
        priceMC = optionMC_E[idx] / P0Td(T)
        IVMC[idx] = ImpliedVolatilityBlack76(CP, priceMC, k, T, frwdFX) * 100.0

    fig2=plt.figure(2)
    plt.plot(K, IVMC)
    plt.plot(K, IVCos, '--r')

    plt.grid()
    plt.legend(['IV-COS', 'IV-MC'])
    plt.title("Fx Implied volatilities")
    plt.xlabel("Strike")
    plt.ylabel("Implied volatility")
    st.pyplot(fig2)




if __name__ == "__main__":
        mainCalculation()