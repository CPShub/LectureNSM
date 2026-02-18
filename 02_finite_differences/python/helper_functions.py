import matplotlib.pyplot as plt
import numpy as np

def plot_force(x:np.array, phi:np.array, title:str, ytitle:str) -> None:
    fig, ax1 = plt.subplots(ncols=1, figsize=[7, 2.5])
    fig.suptitle(title)

    c = "#C24C4C"
    ax1.set_ylim([0, 1.5])

    if "Nm" in ytitle:
        c = "#435384"
        ax1.set_ylim([-1, 0.5])

    ax1.plot(x, phi, label="", color=c, linestyle = "-", marker="o")
    ax1.fill_between(x, phi, 0, label="", color=c, alpha=0.3)
    ax1.set_ylabel(ytitle)
    ax1.set_xlabel("x in m")
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

def plot_deflection(x:np.array, w_numeric:np.array, w_analytic:np.array, title:str, ytitle:str) -> None:
    fig, ax1 = plt.subplots(ncols=1, figsize=[7, 2.5])
    fig.suptitle(title)

    ax1.plot(x, w_analytic, label="Analytische Lösung", color="#435384", linestyle = "-", marker="o")
    ax1.plot(x, w_numeric, label="Numerische Lösung", color="#C24C4C", linestyle = "-", marker="d")
    ax1.set_ylabel(ytitle)
    ax1.set_xlabel("x in m")
    ax1.yaxis.set_inverted(True)
    ax1.legend(loc="upper right")
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

def plot_convergence(step_sizes:np.array, err_a:np.array, err_b:np.array, label_a:str, label_b:str) -> None:
    fig, (ax) = plt.subplots(ncols=1, figsize=[5, 4.5])
    fig.suptitle("Konvergenzanalyse der FD-Methoden")

    ax.loglog(step_sizes, err_a[:, 0], label=label_a, color="#435384", marker='o', linestyle='-')
    ax.loglog(step_sizes, err_b[:, 0], label=label_b, color="#C24C4C", marker='d', linestyle='-')
    ax.set_xlabel("Zeitschrittweite")
    ax.set_ylabel("Abweichung vom analytischen Wert")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(visible=True, which='both', linestyle='--', linewidth=0.5)
    ax.legend(loc = "upper left")

def plot_heat(heat_source:np.array, X:np.array, Y:np.array, Phi:np.array) -> None:
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=[11, 4])
    fig = plt.figure(figsize =(4.5, 4.5))

    heat = ax1.pcolormesh(X, Y, heat_source, cmap = plt.get_cmap('copper'), shading='nearest')
    ax1.axis('equal')
    ax1.set_title("Spezifischer Wärmestrom in W/kg")
    cbar = fig.colorbar(heat, ax=ax1)
    ax1.set_xlabel("x in m")
    ax1.set_ylabel("y in m")

    temp = ax2.pcolormesh(X, Y, Phi, cmap = plt.get_cmap('hot'), shading='nearest')
    cbar = fig.colorbar(temp, ax=ax2)
    ax2.axis('equal')
    ax2.set_title("Temperatur in K")
    ax2.set_xlabel("x in m")
    ax2.set_ylabel("y in m")