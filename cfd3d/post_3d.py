"""Niveau 3 : post-traitement des calculs Code_Saturne (cfd3d/resultats/) et comparaison au 0D.

Usage : python cfd3d/post_3d.py   (après cfd3d/campagne.sh ; écrit figures/3d_*.png)
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RACINE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RACINE))
sys.stdout.reconfigure(encoding="utf-8")
from piscine.modele0d import Piscine, simuler  # noqa: E402
from piscine.physique import G, P_ATM, T_sat, rho_eau  # noqa: E402

RES, FIG = RACINE / "cfd3d" / "resultats", RACINE / "figures"
CAS = ["uniforme", "bord", "coin"]
RACKS = dict(x=(1.0, 11.0), y=(1.0, 9.0), z=(0.3, 4.5))
H = 12.5


def series(nom):
    return np.genfromtxt(RES / nom / "piscine_series.csv", delimiter=",", names=True)


def champ(nom, t=None):
    """Assemble les fichiers par rang du dernier instantané (ou de l'instant t, s) en grilles 3D."""
    fichiers = sorted((RES / nom).glob("champ_t*_r*.csv"))
    instants = sorted({f.name.split("_r")[0] for f in fichiers})
    cle = instants[-1] if t is None else min(instants, key=lambda k: abs(float(k[7:-1]) - t))
    d = np.vstack([np.loadtxt(f, delimiter=",", skiprows=1) for f in fichiers if f.name.startswith(cle + "_")])
    axes = [np.unique(d[:, k]) for k in range(3)]
    idx = [np.searchsorted(axes[k], d[:, k]) for k in range(3)]
    grilles = {}
    for k, n in enumerate(("T", "u", "v", "w")):
        g = np.full([a.size for a in axes], np.nan)
        g[tuple(idx)] = d[:, 3 + k]
        grilles[n] = g
    return float(cle[7:-1]), axes, grilles


# --- 1. Bilan d'énergie et comparaison au 0D ---------------------------------------------------
ref = Piscine(v_solide=0.0)
r0 = simuler(ref, t_max=7 * 3600.0)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4))
a1.plot(r0["t"] / 3600, r0["T"], "k-", lw=2, label="0D")
print("Cas      | écart max T_moy 3D-0D | bilan d'énergie | surchauffe fin | T_surf - T_moy | v_max")
for nom in CAS:
    if not (RES / f"{nom}_nx48").exists():
        continue
    s = series(f"{nom}_nx48")
    t = s["t_s"]
    ecart = np.abs(s["T_moy_C"] - np.interp(t, r0["t"], r0["T"])).max()
    # E(t) - E(0) doit égaler ∫ (P - pertes) dt
    apport = np.concatenate([[0], np.cumsum(np.diff(t) * 0.5 * ((s["P_W"] - s["pertes_W"])[1:]
                                                               + (s["P_W"] - s["pertes_W"])[:-1]))])
    bilan = (s["energie_J"][-1] - s["energie_J"][0]) / apport[-1] - 1
    dT = s["T_max_C"] - s["T_moy_C"]
    print(f"{nom:<8} | {ecart:8.3f} K            | {100*bilan:+6.2f} %        | {dT[-1]:6.2f} K       | "
          f"{s['T_surf_C'][-1]-s['T_moy_C'][-1]:+6.2f} K       | {100*s['v_max_m_s'][-1]:5.1f} cm/s")
    l, = a1.plot(t / 3600, s["T_moy_C"], "--", label=f"3D {nom} : moyenne")
    a1.plot(t / 3600, s["T_max_C"], "-", lw=0.8, c=l.get_color(), label=f"3D {nom} : max")
    a2.plot(t / 3600, dT, c=l.get_color(), label=nom)
a1.set(xlabel="temps (h)", ylabel="T (°C)", title="Code_Saturne 3D vs modèle 0D")
a2.set(xlabel="temps (h)", ylabel="T_max − T_moy (K)", title="Point chaud selon la répartition des assemblages")
for a in (a1, a2):
    a.grid(alpha=0.3)
    a.legend(fontsize=7)
fig.tight_layout()
fig.savefig(FIG / "3d_vs_0d.png", dpi=150)

# --- 2. Coupes : plan vertical médian et plan horizontal au sommet des racks --------------------
cas_dispo = [c for c in CAS if (RES / f"{c}_nx48").exists()]
fig, axs = plt.subplots(2, len(cas_dispo), figsize=(5.2 * len(cas_dispo), 9), squeeze=False)
print("\nPoint chaud (dernier instantané) :")
for j, nom in enumerate(cas_dispo):
    t, (x, y, z), g = champ(f"{nom}_nx48")
    T = g["T"]
    i, k, m = np.unravel_index(np.nanargmax(T), T.shape)
    Tsat = 100.0 + T_sat(P_ATM + rho_eau(30.0) * G * (H - z[m])) - T_sat(P_ATM)
    print(f"  {nom:<8} T_max = {T[i, k, m]:.2f} °C en (x, y, z) = ({x[i]:.1f}, {y[k]:.1f}, {z[m]:.1f}) m ; "
          f"saturation locale {Tsat:.0f} °C, T_moy = {np.nanmean(T):.2f} °C")
    # plan vertical passant par le point chaud
    a = axs[0, j]
    im = a.pcolormesh(x, z, T[:, k, :].T, cmap="inferno", shading="nearest")
    s = 3
    a.quiver(x[::s], z[::s], g["u"][::s, k, ::s].T, g["w"][::s, k, ::s].T, color="w", scale=1.5, width=0.004)
    a.add_patch(plt.Rectangle((RACKS["x"][0], RACKS["z"][0]), 10, RACKS["z"][1] - RACKS["z"][0],
                              fill=False, ec="c", ls="--"))
    a.set(title=f"{nom} – plan y = {y[k]:.1f} m, t = {t/3600:.1f} h", xlabel="x (m)", ylabel="z (m)", aspect="equal")
    fig.colorbar(im, ax=a, shrink=0.8, label="T (°C)")
    # plan horizontal juste au-dessus des racks
    m2 = np.searchsorted(z, RACKS["z"][1])
    a = axs[1, j]
    im = a.pcolormesh(x, y, T[:, :, m2].T, cmap="inferno", shading="nearest")
    a.quiver(x[::s], y[::s], g["u"][::s, ::s, m2].T, g["v"][::s, ::s, m2].T, color="w", scale=1.0, width=0.004)
    a.add_patch(plt.Rectangle((RACKS["x"][0], RACKS["y"][0]), 10, 8, fill=False, ec="c", ls="--"))
    a.set(title=f"plan z = {z[m2]:.1f} m (sommet des racks)", xlabel="x (m)", ylabel="y (m)", aspect="equal")
    fig.colorbar(im, ax=a, shrink=0.8, label="T (°C)")
fig.suptitle("Code_Saturne 3D, k-ω SST : température et vitesses (racks en pointillés, assemblages chauds ×3)")
fig.tight_layout()
fig.savefig(FIG / "3d_coupes.png", dpi=130)

# --- 3. Convergence de maillage (cas uniforme, t = 2 h) -----------------------------------------
print("\nConvergence de maillage (uniforme, t = 2 h) :")
for nx in (32, 48, 64):
    nom = f"uniforme_nx{nx}"
    if not (RES / nom).exists():
        continue
    s = series(nom)
    a2h = lambda c: np.interp(7200.0, s["t_s"], s[c])
    # surchauffe moyennée sur la 2e heure : le panache turbulent fluctue de ±0,03 K
    m = (s["t_s"] > 3600) & (s["t_s"] <= 7200)
    print(f"  nx = {nx:3d} ({nx*round(nx*10/12)*round(nx*12.5/12):>7d} mailles) : T_moy {a2h('T_moy_C'):.3f} °C, "
          f"surchauffe {(s['T_max_C'] - s['T_moy_C'])[m].mean():.3f} K, v_max {100*s['v_max_m_s'][m].mean():.2f} cm/s")
print("\nFigures :", FIG / "3d_vs_0d.png", FIG / "3d_coupes.png")
