import numpy as np
import matplotlib.pyplot as plt

"""
notation we will be using :
Lx,Ly = length and breadth of the 2-dimensional box respectively
n_x,n_y = energy levels
x1,x2 = range of x
y1,y2 = range of y
N = no. of trials
"""

def psi(x,y,n_x,n_y,Lx,Ly):
    # wave function of a particle in a 2D infinite well
    return np.sqrt(4/(Lx*Ly)) * np.sin(n_x * np.pi * x / Lx) * np.sin(n_y * np.pi * y / Ly)

def density(Psi_xy):
    # born probability density |Psi|^2
    return np.abs(Psi_xy)**2

def psi_1d(x, n, L):
    # 1D infinite-well wave function 
    return np.sqrt(2 / L) * np.sin(n * np.pi * x / L)

def mc_probability(n_x, n_y, Lx, Ly, x_range, y_range, N):
    # Monte Carlo estimate of P(x1<=x<=x2, y1<=y<=y2).
    # Returns (estimate, standard_error).
    (x1,x2) = x_range
    (y1,y2) = y_range
    x = np.random.uniform(x1, x2, size=N)
    y = np.random.uniform(y1, y2, size=N)
    psi_xy = psi(x,y,n_x,n_y, Lx, Ly)
    prob_density = density(psi_xy)
    avg_density = prob_density.mean()
    area = abs(x2 - x1) * abs(y2 - y1)
    # by definition, avg_density = int(prob_density dxdy)/[int(dx)*int(dy)] = P/[(x2-x1)*(y2-y1)]
    # so, P = avg_density*(x2-x1)*(y2-y1) = avg_density * area
    P = avg_density * area
    # for error, est. error = sigma/sqrt(N) where sigma = standard deviation
    density_dev = prob_density.std(ddof=1)
    sigma = area * density_dev
    error = sigma/np.sqrt(N)
    return P, error

# we know X and Y are independent variables, since |Psi(x,y)|^2 = |g(x)|^2*|h(y)|^2
# i.e f_X,Y(X=x,Y=y) = f_X(X=x)*f_Y(Y=y) for all x and y [ f is the probability distribution function]
# so P(x1<=x<=x2,y1<=y<=y2) = P(x1<=x<=x2)*P(y1<=y<=y2)
def mc_probability_factorized(x_range, y_range, N, n_x, n_y, Lx, Ly):
    def mc_1d(a, b, n, L):
        t = np.random.uniform(a,b, size = N)
        return (b - a) * np.mean(psi_1d(t, n, L) ** 2)
    P_x = mc_1d(*x_range, n_x, Lx)
    P_y = mc_1d(*y_range, n_y, Ly)
    P = P_x*P_y
    return P

def exact_1d(a, b, n, L):
    # exact integral of |psi_n(t)|^2 from a to b
    # antiderivative: F(t) = t/L - sin(2*n*pi*t/L) / (2*n*pi)
    F = lambda t: t / L - np.sin(2 * n * np.pi * t / L) / (2 * n * np.pi)
    return F(b) - F(a)

def exact_probability(x_range, y_range, n_x, n_y, Lx, Ly):
    # x and y are independent, so the exact probability factorizes too
    return exact_1d(*x_range, n_x, Lx) * exact_1d(*y_range, n_y, Ly)

def convergence(sizes, trials, n_x, n_y, Lx, Ly, x_range, y_range):
    # mean absolute error of the monte carlo estimate for each sample size in sizes
    exact = exact_probability(x_range, y_range, n_x, n_y, Lx, Ly)
    mean_errors = []
    for N in sizes:
        errs = []
        for a in range(trials):
            P, a = mc_probability(n_x, n_y, Lx, Ly, x_range, y_range, N)
            errs.append(abs(P - exact))
        mean_errors.append(np.mean(errs))
    return np.array(mean_errors)

def plot_heatmap(n_x, n_y, Lx, Ly, x_range, y_range, exact):
    # plotting heatmap of |Psi|^2 over the 2d box, with the integration region outlined in red
    g = np.linspace(0, Lx, 300)
    h = np.linspace(0, Ly, 300)
    X, Y = np.meshgrid(g, h)
    D = density(psi(X, Y, n_x, n_y, Lx, Ly))

    fig, ax = plt.subplots()
    im = ax.imshow(D, origin="lower", extent=[0, Lx, 0, Ly], cmap="viridis")
    fig.colorbar(im, ax=ax, label="|Psi|^2")

    (x1, x2), (y1, y2) = x_range, y_range
    ax.add_patch(plt.Rectangle((x1, y1), x2 - x1, y2 - y1,
                               fill=False, edgecolor="red", linewidth=2,
                               label=f"region, P = {exact:.4f}"))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Probability density")
    ax.legend(loc="upper right")
    fig.savefig("heatmap.png", dpi=150)
    plt.show()

if __name__ == "__main__":
    # we would like our results to be reproducible
    # hence picking a seed = 42 (Answer to the Ultimate Question of Life, the Universe, and Everything)
    np.random.seed(42)
    Lx, Ly = 1, 1
    n_x, n_y = 1, 1
    x_range = (0, Lx / 2)
    y_range = (0, Ly / 2)
    N = 1_000_000

    # checking w.r.t exact formula
    assert np.isclose(exact_1d(0, Lx, n_x, Lx), 1)      # normalization
    assert np.isclose(exact_1d(0, Lx / 2, 1, Lx), 0.5)  # symmetry, n = 1

    exact = exact_probability(x_range, y_range, n_x, n_y, Lx, Ly)
    plot_heatmap(n_x, n_y, Lx, Ly, x_range, y_range, exact)
    P, error = mc_probability(n_x, n_y, Lx, Ly, x_range, y_range, N)
    P_fact = mc_probability_factorized(x_range, y_range, N, n_x, n_y, Lx, Ly)

    print(f"Exact              : {exact:.6f}")
    print(f"Monte Carlo (2D)            : {P:.6f} +/- {error:.6f}")
    print(f"Monte Carlo (P_x * P_y)     : {P_fact:.6f}")

    # use 9 sample sizes spaced evenly on a log scale, from 10^2 to 10^6
    sizes = np.logspace(2, 6, 9).astype(int)
    # each N is run 50 times and the absolute errors are averaged since a single iteration is too noisy to show a trend
    errors = convergence(sizes, 50, n_x, n_y, Lx, Ly, x_range, y_range)

    plt.loglog(sizes, errors, "o-", label="Monte Carlo error (mean of 50 trials)")
    plt.loglog(sizes, errors[0] * np.sqrt(sizes[0] / sizes), "k--", label="1/sqrt(N) reference")
    plt.xlabel("N (samples)")
    plt.ylabel("absolute error = |estimate - exact|")
    plt.legend()
    plt.savefig("mc_error.png", dpi=150)
    plt.show()
