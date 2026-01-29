
import matplotlib.pyplot as plt
import numpy as np

def plot_motion(t:np.array, x_lin:np.array, x_nlin:np.array, title:str) -> None:
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=[14, 2.5])
    fig.suptitle(title)

    ax1.plot(t, x_lin[:, 0], label="Linearisiert", color="#008374", linestyle="-")
    ax1.plot(t, x_nlin[:, 0], label="Nicht Linearisiert", color="#344b47", linestyle="-")
    ax1.set_ylabel("Φ in rad")
    ax1.set_xlabel("t in s")
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    ax2.plot(t, x_lin[:, 1], label="Linearisiert", color="#008374", linestyle="-")
    ax2.plot(t, x_nlin[:, 1], label="Nicht Linearisiert", color="#344b47", linestyle="-")
    ax2.set_ylabel("v in rad/s")
    ax2.set_xlabel("t in s")
    ax2.legend(loc = "upper right")
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

def plot_energy(t:np.array, E_pot:list, E_kin:list, E_total:list, E_pot_nonlinear:list, E_kin_nonlinear:list, E_total_nonlinear:list, title:str) -> None:
    fig4, (ax7, ax8) = plt.subplots(ncols=2, figsize=[14, 2.5])
    fig4.suptitle(title)

    ax7.plot(t, E_pot, label="Pot. Energie", color="#008374", linestyle="-")
    ax7.plot(t, E_kin, label="Kin. Energie", color="#936496", linestyle="-")
    ax7.plot(t, E_total, label="Ges. Energie", color="#344b47", linestyle="-")
    ax7.set_ylabel("Energie in J")
    ax7.set_xlabel("t in s")
    ax7.spines['top'].set_visible(False)
    ax7.spines['right'].set_visible(False)
    ax7.set_title("Linearisiert", y=0.93)


    ax8.plot(t, E_pot_nonlinear, label="Pot. Energie", color="#008374", linestyle="-")
    ax8.plot(t, E_kin_nonlinear, label="Kin. Energie", color="#936496", linestyle="-")
    ax8.plot(t, E_total_nonlinear, label="Ges. Energie", color="#344b47", linestyle="-")
    ax8.set_ylabel("Energie in J")
    ax8.set_xlabel("t in s")
    ax8.set_title("Nicht Linearisiert", y=0.93)
    ax8.spines['top'].set_visible(False)
    ax8.spines['right'].set_visible(False)
    ax8.legend()

def plot_convergence(timesteps:list, error_se:list, error_tr:list) -> None:
    fig, ax = plt.subplots(ncols=1, figsize=[4.5, 4.5])

    ax.loglog(timesteps, error_se, label='Symplektisches Euler-Verfahren', color="#008374")
    ax.loglog(timesteps, error_tr, label='Trapezregel', color="#344b47")
    ax.set_xlabel("Zeitschrittgröße Δt in s")
    ax.set_ylabel("Abweichung vom analytischen Wert in rad")
    ax.set_title("Konvergenzanalyse der Zeitintegrationsverfahren")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(visible=True, which='both', linestyle='--', linewidth=0.5)
    ax.legend()