import matplotlib.pyplot as plt
import numpy as np

def plot_heat(heat_source:np.array, X:np.array, Y:np.array, Phi:np.array) -> None:
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=[11, 4])

    heat = ax1.pcolormesh(X, Y, heat_source, cmap = plt.get_cmap('copper'), shading='nearest')
    ax1.axis('equal')
    ax1.set_title("Spezifischer Wärmestrom in W/m²")
    cbar = fig.colorbar(heat, ax=ax1)
    ax1.set_xlabel("x in m")
    ax1.set_ylabel("y in m")

    temp = ax2.pcolormesh(X, Y, Phi, cmap = plt.get_cmap('hot'), shading='nearest')
    cbar = fig.colorbar(temp, ax=ax2)
    ax2.axis('equal')
    ax2.set_title("Temperatur in K")
    ax2.set_xlabel("x in m")
    ax2.set_ylabel("y in m")