
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from IPython.display import HTML

def plot_motion(t:np.array, x_ana:np.array, x_lin:np.array, x_nlin:np.array, title:str) -> None:
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=[14, 2.5])
    fig.suptitle(title)

    ax1.plot(t, x_ana[:, 0], label="Analytisch", color="#435384", linestyle="--")
    ax1.plot(t, x_lin[:, 0], label="Linearisiert", color="#435384", linestyle="-")
    ax1.plot(t, x_nlin[:, 0], label="Nicht Linearisiert", color="#C24C4C", linestyle="-")
    ax1.set_ylabel("Φ in rad")
    ax1.set_xlabel("t in s")
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    ax2.plot(t, x_ana[:, 1], label="Analytisch", color="#435384", linestyle="--")
    ax2.plot(t, x_lin[:, 1], label="Linearisiert", color="#435384", linestyle="-")
    ax2.plot(t, x_nlin[:, 1], label="Nicht Linearisiert", color="#C24C4C", linestyle="-")
    ax2.set_ylabel("v in rad/s")
    ax2.set_xlabel("t in s")
    ax2.legend(loc = "upper right")
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

def plot_energy(t:np.array, E_pot:list, E_kin:list, E_total:list, E_pot_nonlinear:list, E_kin_nonlinear:list, E_total_nonlinear:list, title:str) -> None:
    fig4, (ax7, ax8) = plt.subplots(ncols=2, figsize=[14, 2.5])
    fig4.suptitle(title)

    ax7.plot(t, E_pot, label="Pot. Energie", color="#435384", linestyle="-")
    ax7.plot(t, E_kin, label="Kin. Energie", color="#C24C4C", linestyle="-")
    ax7.plot(t, E_total, label="Ges. Energie", color="#F6A315", linestyle="-")
    ax7.set_ylabel("Energie in J")
    ax7.set_xlabel("t in s")
    ax7.spines['top'].set_visible(False)
    ax7.spines['right'].set_visible(False)
    ax7.set_title("Linearisiert", y=0.93)


    ax8.plot(t, E_pot_nonlinear, label="Pot. Energie", color="#435384", linestyle="-")
    ax8.plot(t, E_kin_nonlinear, label="Kin. Energie", color="#C24C4C", linestyle="-")
    ax8.plot(t, E_total_nonlinear, label="Ges. Energie", color="#F6A315", linestyle="-")
    ax8.set_ylabel("Energie in J")
    ax8.set_xlabel("t in s")
    ax8.set_title("Nicht Linearisiert", y=0.93)
    ax8.spines['top'].set_visible(False)
    ax8.spines['right'].set_visible(False)
    ax8.legend()

def plot_convergence(timesteps:list, error_se:list, error_mp:list, error_tr:list) -> None:
    fig, ax = plt.subplots(figsize=[7, 5])
    fig.suptitle("Konvergenzanalyse der Zeitintegrationsverfahren")

    ax.loglog(timesteps, error_se, label='Symplektisches Euler-Verfahren', color="#435384", marker='o', linestyle='-')
    ax.loglog(timesteps, error_mp, label='Explizite Mittelpunktsregel', color="#F6A315", marker='o', linestyle='-')
    ax.loglog(timesteps, error_tr, label='Trapezregel', color="#C24C4C", marker='o', linestyle='-')
    ax.set_xlabel("Zeitschrittweite Δt in s")
    ax.set_ylabel("E = |Phi - Phi_ref| [rad]")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(visible=True, which='both', linestyle='--', linewidth=0.5)
    ax.legend(loc = "upper left")

def make_animation(t:np.array, vals_lin:np.array, vals_nlin:np.array, bounds:float, delta_t:float, title:str):
    fig, (ax, ax2) = plt.subplots(ncols=2, figsize=[10, 4])
    fig.suptitle(title)
    ax.autoscale(False)
    ax.set_xlim((-bounds, bounds))
    ax.set_ylim((-bounds, bounds))
    ax.set_aspect('equal')
    ax.set_title("Linearisiert")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    line, = ax.plot([], [], 'o-', lw=2, color='#435384')

    ax2.autoscale(False)
    ax2.set_xlim((-bounds, bounds))
    ax2.set_ylim((-bounds, bounds))
    ax2.set_aspect('equal')
    ax2.set_title("Nicht Linearisiert")
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    line2, = ax2.plot([], [], 'o-', lw=2, color='#C24C4C')
    time_template = '%.1f s'
    time_text = ax2.text(0.05, 0.9, '', transform=ax2.transAxes)

    def animate(i):
        x_lin = np.sin(vals_lin[i, 0])
        y_lin = -np.cos(vals_lin[i, 0])
        x_nlin = np.sin(vals_nlin[i, 0])
        y_nlin = -np.cos(vals_nlin[i, 0])
        line.set_data([0, x_lin], [0, y_lin])
        time_text.set_text(time_template % (i*delta_t))
        line2.set_data([0, x_nlin], [0, y_nlin])
        return line, line2, time_text

    ani = animation.FuncAnimation(
        fig, animate, frames=len(t), interval=delta_t*1000, blit=True
    )
    plt.close()

    return HTML(ani.to_jshtml(fps=20))